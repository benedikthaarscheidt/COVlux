
%% 0) Gurobi threads & paths
nThreads = feature('numcores')-2; %dont use all cores since computer might crash
fprintf('Using %d threads for Gurobi\n',nThreads);
setenv('OMP_NUM_THREADS',num2str(nThreads));
addpath(fullfile(getenv('GUROBI_HOME'),'examples','matlab'),'-begin');

%% 1) Load irreversible pruned model
data     = load('/Users/benedikthaarscheidt/M.Sc./master_thesis_second_moment/Models/generic_models/E_coli/iML1515_pruned_permissive_biomass_loopless_v6.mat','pruned_ir')%iML1515_pruned_permissive_biomassrm_loopless_v5
model_ir = data.pruned_ir;
S        = model_ir.S;
rxnNames = model_ir.rxns;
[m,n]    = size(S);

met_indices = find(model_ir.S(:, 2219) ~= 0);

% 2. Extract the coefficients and the metabolite names/IDs
coefficients = full(model_ir.S(met_indices, 2219));
met_ids = model_ir.mets(met_indices);

% 3. Try to get the descriptive names if they exist in the model
if isfield(model_ir, 'metNames')
    met_names = model_ir.metNames(met_indices);
else
    met_names = met_ids; % Fallback to IDs
end

% 4. Display as a table for easy reading
biomass_structure = table(met_ids, met_names, coefficients, ...
    'VariableNames', {'ID', 'Name', 'Stoichiometry'});

disp(biomass_structure);

%% 2) Parameters
eps_flux      = 1e-7;
tol_balance   = 1e-8;
M             = 5e3;
badPair_count = 0;
JACCARD_TAU   = 0.99;              % skip a candidate if Jaccard(similarity) >= 0.95 to any accepted support
MAX_PER_ANCHOR = 40; % accept up to this many supports per anchor (set Inf for no cap)
MIN_PER_ANCHOR=3;

accepted_supps = {};    % global list of supports to compare Jaccard against
EFM_supps  = {};        % all accepted supports (any anchor)
EFM_vecs   = {};        % corresponding full v vectors
EFM_anchor = [];        % anchor index for each accepted support

%% Compute range of all reactions for fixing      the flux 

fvaModel.A        = sparse(S);
fvaModel.rhs      = zeros(m,1);
fvaModel.sense    = repmat('=', m, 1);
fvaModel.lb       = model_ir.lb;
fvaModel.ub       = model_ir.ub;
fvaModel.vtype    = repmat('C', n, 1);
fvaModel.modelsense = 'min';%will be switched to max
fvaModel.obj      = zeros(n,1);

% Gurobi parameters for FVA
fvaParams.OutputFlag = 0;
fvaParams.TimeLimit  = 60;    
fvaParams.FeasibilityTol = 1e-9;

vMin = nan(n,1);
vMax = nan(n,1);

for j = 1:n
    % 1) minimize v_j
    fvaModel.obj(:) = 0;
    fvaModel.obj(j) = 1;
    fvaModel.modelsense = 'min';
    sol = gurobi(fvaModel, fvaParams);
    if isfield(sol,'status') && strcmp(sol.status,'OPTIMAL')
        vMin(j) = sol.objval;
    else
        warning('FVA min infeasible for reaction %d (%s)', j, rxnNames{j});
    end

    % 2) maximize v_j
    fvaModel.modelsense = 'max';
    sol = gurobi(fvaModel, fvaParams);
    if isfield(sol,'status') && strcmp(sol.status,'OPTIMAL')
        vMax(j) = sol.objval;
    else
        warning('FVA max infeasible for reaction %d (%s)', j, rxnNames{j});
    end
end

range=[vMin,vMax];


%%
tol0 = 1e-5;                        
idxZeroMax = abs(range(:,2)) <= tol0;

fprintf('Reactions with vMax <= %.1e: %d\n', tol0, nnz(idxZeroMax));
disp(rxnNames(idxZeroMax));

if any(idxZeroMax)
    T = table(rxnNames(idxZeroMax), range(idxZeroMax,1), range(idxZeroMax,2), ...
              'VariableNames', {'rxn','vMin','vMax'});
    disp(T);
end
%% 3) Build static MILP template 

vMin(~isfinite(vMin)) = 0;          % irreversibles: safe floor
vMax(~isfinite(vMax)) = M;          % cap unbounded
vMin = max(vMin, 0);                % safety for irreversible split model
range = [vMin, vMax];

lb = [range(:,1); zeros(n,1)];
ub = [range(:,2);  ones(n,1)];
ub_eff = ub; 
ub_eff(~isfinite(ub_eff)) = M;
ub = ub_eff;
%%
vtype = [ repmat('C',n,1) ; repmat('B',n,1) ];


Aeq = [ sparse(S), sparse(m,n) ];
beq = zeros(m,1);

%    v_i - M*b_i <= 0
%   -v_i + eps*b_i <= 0
A1 = [  speye(n),  -ub(1:n).*speye(n)
       -speye(n),  eps_flux*speye(n) ];
b1 = zeros(2*n,1);


nz         = (S~=0);
importCols = false(1,n);
exportCols = false(1,n);
for j = 1:n
    rows = nz(:,j);
    if any(rows)
        col = S(rows,j);
        importCols(j) = all(col> 0);
        exportCols(j) = all(col< 0);
    end
end
Aimp = []; bimp = [];
if any(importCols)
    r = sparse(1,2*n);  
    r(n+find(importCols)) = -1;
    Aimp = r;
    bimp = -1;
end


Aexp = []; bexp = [];
if any(exportCols)
    r = sparse(1,2*n);
    r(n+find(exportCols)) = -1;
    Aexp = r;
    bexp = -1;
end




%    b_i + b_j <= 1
badPartner = zeros(n,1);
badPair_count = 0;
for i = 1:n
    rn = rxnNames{i};
    if endsWith(rn,'_f')
        j = find(strcmp(rxnNames,[rn(1:end-2) '_b']),1);
        if ~isempty(j)
            badPartner(i)=j;
            badPartner(j)=i;
            if i < j  % count each pair only once
                badPair_count = badPair_count + 1;
            end
        end
    end
end
disp(badPair_count)

bpairs = find(badPartner>0 & (1:n)' < badPartner);
p = numel(bpairs);
Abp = sparse(p,n*2);
bbp = ones(p,1);
for k=1:p
    i = bpairs(k);
    j = badPartner(i);
    Abp(k, n+i) = 1;
    Abp(k, n+j) = 1;
end

fprintf('IMPORT reactions classified (%d):\n', sum(importCols));
fprintf('Export reactions classified (%d):\n', sum(exportCols));

Aineq = [ A1
          Aimp
          Aexp
          Abp];
bineq = [ b1
          bimp
          bexp 
          bbp];

sense = [ repmat('=', size(Aeq,1),1)
          repmat('<', size(Aineq,1),1) ];

obj = [ zeros(n,1) ; ones(n,1) ];


linModel.A          = [Aeq; Aineq];
linModel.rhs        = [beq; bineq];
linModel.sense      = sense;
linModel.lb         = lb;
linModel.ub         = ub;
linModel.vtype      = vtype;
linModel.modelsense = 'min';
linModel.obj        = obj;

%% 4) Gurobi parameters
commonParams = struct( ...
   'FeasibilityTol',1e-9, ...
   'IntFeasTol',    1e-9, ...
   'NumericFocus',  3);
pPool = struct( ...
  'OutputFlag',     0, ...
  'PoolSearchMode',2, ...
  'PoolSolutions',100, ...
  'TimeLimit',    60, ...
  'MIPGap',0.0, ...
  'Cuts',1, ...
  'Presolve',2, ...%this was 2 
  'Heuristics',0.5, ...
  'MIPFocus',1, ...
  'Threads',nThreads );
pSeed = struct( ...g
  'OutputFlag',     0, ...
  'PoolSearchMode',2, ...
  'PoolSolutions',40, ...
  'TimeLimit',    600, ...
  'MIPGap',0.3, ...
  'Cuts',-1, ... 
  'Presolve',2, ...
  'Heuristics',0.2, ...
  'MIPFocus',1, ...
  'Threads',nThreads );
for p = {pPool,pSeed}
    p{1}.FeasibilityTol = commonParams.FeasibilityTol;
    p{1}.IntFeasTol     = commonParams.IntFeasTol;
    p{1}.NumericFocus   = commonParams.NumericFocus;
end

fvaParams.OutputFlag   = 0;    
fvaParams.TimeLimit    = 60;  
fvaParams.FeasibilityTol = 1e-9; 



%% 6) Enumeration loop
covered  = false(n,1);
skipped  = false(n,1);
allEFMs  = {};
allX     = {};
prevUc   = 0;
efm_count_total       = 0;              % how many EFMs accepted
efm_with_badpair_cnt  = 0; 
p = numel(bpairs);% how many EFMs contain >=1 bad pair
badpair_counts        = zeros(p,1);

%% === Enumeration loop: accept-if-fits, else prune-to-3-checks ===

maxSeedAttempts = 25;
stallLimit = 300;
noProgressCount   = 0;
prevCoveredCount  = sum(covered);
stallTerminated   = false;


pool = gcp('nocreate');
if isempty(pool), pool = parpool('local', nThreads); end

%% 5) Checkpoint / Resume (timed saves only)
CHECKPOINT_ENABLE  = true;
AUTOSAVE_SEC       = 15*60;                      % save every 15 minutes
CHECKPOINT_FILE    = fullfile(pwd,'efm_checkpoint.mat');

% Optional: resume from an existing checkpoint file
LOAD_CHECKPOINT = false;                         % set true to resume
LOAD_FILE       = CHECKPOINT_FILE;               % or point to a specific file

% Lightweight model signature to avoid mismatched resumes
model_sig = struct('m',m,'n',n,'nnzS',nnz(S),'rxnCount',numel(rxnNames));

if LOAD_CHECKPOINT && exist(LOAD_FILE,'file')
    R = load(LOAD_FILE,'model_sig');
    if isfield(R,'model_sig') && isequal(R.model_sig, model_sig)
        fprintf('[resume] Loading: %s\n', LOAD_FILE);
        load_checkpoint_state(LOAD_FILE);        % helper at end of file
        if exist('rng_state','var'), rng(rng_state); end
    else
        warning('[resume] Incompatible checkpoint; starting fresh.');
    end
end


if ~exist('covered','var'), covered = false(n,1); end
if ~exist('skipped','var'), skipped = false(n,1); end

% Timer for periodic saves
lastSave = tic;
targetBiomassName = 'BIOMASS_Ec_iML1515_WT_75p37M';
bio_idx = find(strcmp(model_ir.rxns, targetBiomassName));
fprintf('Biomass_idx: %d\n', bio_idx);

%% 6) Enumeration loop (timed checkpoint only; logic unchanged)
try
while any(~covered & ~skipped)
    % Timed checkpoint (even if no acceptance happens in this iteration)
    if CHECKPOINT_ENABLE && toc(lastSave) > AUTOSAVE_SEC
        save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
        lastSave = tic;
    end

    toTry = find(~covered & ~skipped);
    fprintf('\n=== Next reaction: %d remaining (%d covered) ===\n', numel(toTry), sum(covered));
    idx = randi(numel(toTry));
    idx=bio_idx;
    i = toTry(idx);
    fprintf('I: %d\n', i);
    prevUc = i;

    localTemplate = linModel;
    found = false;
    accepted_this_anchor = 0;   

    for attempt = 1:maxSeedAttempts
        M1 = localTemplate;
        M1.lb(n+i) = 1;
        M1.ub(n+i) = 1;
        M1.lb(i)=1;

        sol1 = gurobi(M1, pSeed);
          if ~isfield(sol1,'pool') || isempty(sol1.pool)
            fprintf('  attempt %d: no pool solutions → retrying\n', attempt);
            continue;
        end

        numCand  = numel(sol1.pool);
        fprintf('  attempt %d: %d pool candidates\n', attempt, numCand);

        supp_list        = cell(numCand,1);
        v_list           = cell(numCand,1);
        hasPair_list     = false(numCand,1);
        pairIdx_list     = cell(numCand,1);
        badpair_violate  = false(numCand,1);
        acceptType       = zeros(numCand,1);
        supp_accept      = cell(numCand,1);
        v_accept         = cell(numCand,1);
        res_full_acc     = inf(numCand,1);
        res_supp_acc     = inf(numCand,1);
        maxOff_acc       = zeros(numCand,1);
        idxOff_acc       = zeros(numCand,1);
        nogood_supp      = cell(numCand,1);
        anchorInSupp     = false(numCand,1);

        parfor k = 1:numCand
            xk  = getSolutionVector(sol1.pool(k));
            v_k = xk(1:n);
            supp0 = find(v_k >= eps_flux - tol_balance);
            if isempty(supp0) || ~ismember(i, supp0)
                fprintf("empty support \n")
                continue;
            end
            if ~no_bad_pairs(supp0, badPartner)
                [supp0, v_k] = cancel_twins(supp0, v_k, badPartner, S, eps_flux, tol_balance);
                if isempty(supp0) || ~ismember(i, supp0)
                    continue;
                end
            end
            anchorInSupp(k) = true;
            nogood_supp{k}  = supp0;
            [hBP, pIdx] = badpairs_in_support(supp0, badPartner, bpairs);
            hasPair_list(k) = hBP;
            if hBP, pairIdx_list{k} = pIdx; end

            nullOK      = nullity1_ok(S, supp0);
            res_supp    = norm(S(:,supp0) * v_k(supp0), Inf);
            full_suppOK = (res_supp <= 1e-8);
            stillBP     = ~no_bad_pairs(supp0, badPartner);
            badpair_violate(k) = stillBP;

            if nullOK && full_suppOK && ~stillBP
                v_full = v_k;
                res_full = norm(S*v_full, Inf);
                [res_supp_dbg, maxOff, idxOff] = support_debug_metrics(S, v_full, supp0);
                acceptType(k)   = 1;
                supp_accept{k}  = supp0;
                v_accept{k}     = v_full;
                res_full_acc(k) = res_full;
                res_supp_acc(k) = res_supp_dbg;
                maxOff_acc(k)   = maxOff;
                idxOff_acc(k)   = idxOff;
                continue;
            end

            if (~nullOK || ~full_suppOK)
                fprintf("Nullity: %d, Residual: %d \n",nullOK,full_suppOK)
                [ok_fix, supp_fix, v_fix] = prune_nullity_by_single_drop(S, ub, supp0, eps_flux, tol_balance, badPartner, i);
                if ok_fix
                    res_supp_fix = norm(S(:,supp_fix) * v_fix(supp_fix), Inf);
                    full_suppOK_fix = (res_supp_fix <= 1e-8);
                    stillBP_fix = ~no_bad_pairs(supp_fix, badPartner);
                    if full_suppOK_fix && ~stillBP_fix
                        res_full_fix = norm(S*v_fix, Inf);
                        [res_supp_dbg, maxOff, idxOff] = support_debug_metrics(S, v_fix, supp_fix);
                        acceptType(k)   = 2;
                        supp_accept{k}  = supp_fix(:);
                        v_accept{k}     = v_fix(:);
                        res_full_acc(k) = res_full_fix;
                        res_supp_acc(k) = res_supp_dbg;
                        maxOff_acc(k)   = maxOff;
                        idxOff_acc(k)   = idxOff;
                    end
                end
            end
        end

        for k = 1:numCand
            if ~anchorInSupp(k), continue; end
            if hasPair_list(k) && ~isempty(pairIdx_list{k})
                badpair_counts(pairIdx_list{k}) = badpair_counts(pairIdx_list{k}) + 1;
            end

            if badpair_violate(k)
                efm_with_badpair_cnt = efm_with_badpair_cnt + 1;
                continue;
            end

            if acceptType(k) == 1 || acceptType(k) == 2
                if acceptType(k) == 1
                    supp0 = supp_accept{k};
                    v_k   = v_accept{k};
                    tag   = '';   % normal accept
                else
                    supp0 = supp_accept{k};
                    v_k   = v_accept{k};
                    tag   = ' (single-drop)';
                end

                % --- Jaccard similarity screen (against all previously accepted EFMs) ---
                J = 0;
                if ~isempty(accepted_supps)
                    J = jaccard_max_sim(supp0, accepted_supps);
                end
                if J >= JACCARD_TAU
                    % forbid exact re-selection of this support later
                    localTemplate = applyNoGoodCut(localTemplate, supp0, n);
                    continue;
                end

                % --- Refit: maximize total flux sum on the fixed support ---
                [okLP, v_full2, objval] = maximize_sum_v_on_support(S, ub, supp0, eps_flux);
                if okLP
                    v_full = v_full2;
                else
                    % fall back to candidate’s v on support
                    v_full = zeros(n,1);
                    v_full(supp0) = v_k(supp0);
                end

                % --- Record globally ---
                EFM_supps{end+1}   = supp0(:);
                EFM_vecs{end+1}    = v_full(:);
                EFM_anchor(end+1)  = i;
                accepted_supps{end+1} = supp0(:);

                % coverage + forbid exact support
                covered(supp0)   = true;
                localTemplate    = applyNoGoodCut(localTemplate, supp0, n);

            
                fprintf('   → accept%s |supp|=%d | Sv=%.2e | sum v=%.3g | J=%.2f\n', ...
                        tag, numel(supp0), norm(S*v_full,Inf), sum(v_full(supp0)), J);

                accepted_this_anchor = accepted_this_anchor + 1;
                localTemplate = applyNoGoodCut(localTemplate, nogood_supp{k}, n);

                % --- timed checkpoint after acceptance 
                if CHECKPOINT_ENABLE && toc(lastSave) > AUTOSAVE_SEC
                    save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
                    lastSave = tic;
                end

                
            else
                % reject → add no-good cut so solver moves on
                localTemplate = applyNoGoodCut(localTemplate, nogood_supp{k}, n);
            end
        end % for k

        
        if accepted_this_anchor >= MAX_PER_ANCHOR || accepted_this_anchor>= MIN_PER_ANCHOR
            break;
        end
        
    end % for attempt

    % If nothing accepted for this anchor, mark as skipped 
    if accepted_this_anchor == 0
        warning('   No valid support for reaction %d (%s) after %d attempts', i, rxnNames{i}, maxSeedAttempts);
        skipped(i) = true;
    end

    % Progress / stall logic (unchanged)
    currCoveredCount = sum(covered);
    if currCoveredCount > prevCoveredCount
        noProgressCount = 0;
    else
        noProgressCount = noProgressCount + 1;
    end
    prevCoveredCount = currCoveredCount;

    if noProgressCount >= stallLimit
        fprintf('\nNo new coverage for %d consecutive iterations. Early termination.\n', stallLimit);
        stallTerminated = true;

        % save a final checkpoint before breaking
        if CHECKPOINT_ENABLE
            save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
            lastSave = tic;
        end
        break;
    end
    break
end % while
catch ME
    fprintf('[error] %s\n', ME.message);
    if CHECKPOINT_ENABLE
        try
            save_checkpoint_atomic(CHECKPOINT_FILE, model_sig)
            catch
        end
    end
    rethrow(ME);
end



%%
missingReactions = find(~covered);
fprintf('Uncovered reactions (%d):\n', numel(missingReactions));
disp(rxnNames(missingReactions));

%% 7) Report & save — columns = EFMs, rows = [residual; reaction fluxes]
numEFMs = numel(EFM_vecs);
if numEFMs == 0
    warning('No EFMs accepted; nothing to save.');
else
    % Build matrix: first row = ||S*v||_inf, then n reaction rows
    EFM_matrix = zeros(n+1, numEFMs);
    for c = 1:numEFMs
        v = EFM_vecs{c}(:);              % full-length n×1, zeros off-support
        EFM_matrix(2:end, c) = v;
        EFM_matrix(1, c)     = norm(S*v, Inf);
    end

    % Names
    rowNames = [{'Sv_inf_residual'}; rxnNames(:)];
    rawColNames = arrayfun(@(c) sprintf('EFM_%d_anchor_%s', c, rxnNames{EFM_anchor(c)}), ...
                           1:numEFMs, 'UniformOutput', false);
    varNames = matlab.lang.makeValidName(rawColNames, 'ReplacementStyle','delete');
    varNames = matlab.lang.makeUniqueStrings(varNames);

    % Table view
    EFM_table = array2table(EFM_matrix, 'RowNames', rowNames, 'VariableNames', varNames);

    % Binary support mask (threshold like elsewhere)
    tol_nz = max(eps_flux - tol_balance, 0);
    EFM_support = abs(EFM_matrix(2:end, :)) > tol_nz;

    % Save MAT + CSV (include new structures)
    save('efms_matrix_iML1515_forBiomass2.mat', ...
         'EFM_matrix','EFM_table','rowNames','varNames','EFM_support', ...
         'EFM_anchor','EFM_supps','rxnNames','S','model_ir');

    writetable(EFM_table, ...
               'efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5_higher_cover_v2.csv', ...
               'WriteRowNames', true);

    fprintf('Saved EFM matrix: (%d x %d) = [residual; %d reactions] x %d EFMs\n', ...
            size(EFM_matrix,1), size(EFM_matrix,2), n, numEFMs);
end


%%
load('efms_matrix_iML1515_forBiomass.mat','EFM_matrix','EFM_table','rowNames','varNames','EFM_support','EFM_anchor','EFM_supps','rxnNames','S','model_ir');


if strcmp(EFM_table.Properties.RowNames{1}, 'Sv_inf_residual')
    flux_data = EFM_table{2:end, :};
    rxn_list  = EFM_table.Properties.RowNames(2:end);
else
    flux_data = EFM_table{:, :};
    rxn_list  = EFM_table.Properties.RowNames;
end

num_efms = size(flux_data, 2);
epsilon  = 1e-8; % Threshold for "active" 

% 2. Loop through each EFM and Analyze
fprintf('\nbla=== EFM ANALYSIS REPORT ===\n');
for k = 1:num_efms
    v = flux_data(:, k);
    active_idx = find(v > epsilon);
    support_size = length(active_idx);
    
    fprintf('\n------------------------------------------------------------\n');
    fprintf('EFM #%d (Support Size: %d reactions)\n', k, support_size);
    
    % A. Sanity Check for Biomass EFM
    if support_size < 50
        fprintf('  [WARNING] Support is suspiciously small. Likely a broken/invalid EFM.\n');
    elseif support_size > 200
        fprintf('  [INFO] Support size looks like a valid growth pathway.\n');
    end
    
    % B. Identify Secreted Products (Exports)
    % Logic: Find active Exchange/Sink reactions that CONSUME a metabolite 
    % (meaning they remove it from the system -> Export)
    fprintf('  PRODUCED / EXPORTED MOLECULES:\n');
    found_export = false;
    
    for r_idx = active_idx'
        r_name = rxn_list{r_idx};
        
        % Check if it is an Exchange (EX_), Sink (SK_/DM_), or Biomass
        is_exchange = startsWith(r_name, {'EX_', 'SK_', 'DM_'});
        is_biomass  = contains(r_name, 'BIOMASS');
        
        if is_exchange || is_biomass
            % Check Stoichiometry: 
            % In irreversible models, if an export reaction is active, 
            % it typically has a negative coefficient for the metabolite it removes.
            met_indices = find(S(:, r_idx) < 0); 
            
            if ~isempty(met_indices)
                % Get metabolite names
                if isfield(model_ir, 'metNames')
                    names = model_ir.metNames(met_indices);
                else
                    names = model_ir.mets(met_indices);
                end
                
                % Print the product
                for m = 1:numel(names)
                    fprintf('    -> %s (Flux: %.2e) via %s\n', names{m}, v(r_idx), r_name);
                end
                found_export = true;
            elseif is_biomass
                 fprintf('    -> ** CELL BIOMASS ** (Flux: %.2e)\n', v(r_idx));
                 found_export = true;
            end
        end
    end
    
    if ~found_export
        fprintf('    (No explicit exports found - check Model Definitions)\n');
    end

    % C. List All Active Reactions
    % fprintf('  ACTIVE REACTIONS:\n');
    % disp(rxn_list(active_idx));
end
%% Reconstruct 'covered' from saved EFMs 
% Assumes EFM_matrix (rows: 1 residual + n reactions, cols: EFMs) is loaded.
% Produces: covered  (n x 1 logical), true if reaction active in ≥1 EFM.

% Robust threshold 
if exist('eps_flux','var') && exist('tol_balance','var')
    tol_nz = max(eps_flux - tol_balance, 0);
else
    tol_nz = 1e-6;  % safe default if not in workspace
end

% Handle empty/no EFMs
if ~exist('EFM_matrix','var') || isempty(EFM_matrix)
    n = size(S,2);                % fall back to model size if available
    if isempty(n), error('Cannot infer n without S or EFM_matrix.'); end
    covered = false(n,1);
else
    n = size(EFM_matrix,1) - 1;   % drop residual row
    active = abs(EFM_matrix(2:end, :)) > tol_nz;   % (n x numEFMs) logical
    covered = any(active, 2);                        % (n x 1) logical
end


if exist('rxnNames','var') && numel(rxnNames) ~= numel(covered)
    warning('Size mismatch: rxnNames (%d) vs covered (%d).', numel(rxnNames), numel(covered));
end


%% === EFM basis diagnostics: rank & conditioning ===
% EFM_matrix: (1 + m) x s, first row = ||S*v||_inf, rows 2..m+1 = reaction fluxes
E_full = EFM_matrix(2:end, :);          % m x s (reactions x EFMs)

% remove all-zero columns (shouldn’t happen, but be safe)
col_nz = any(abs(E_full) > 0, 1);
E = E_full(:, col_nz);
s_eff = size(E,2);
m_eff = size(E,1);

if s_eff == 0
    fprintf('\n=== EFM basis diagnostics ===\n');
    fprintf('All columns were zero. Nothing to diagnose.\n');
    return;
end

% Column-normalize to reduce scale effects on numerical rank
colnorm = vecnorm(E, 2, 1);
colnorm(colnorm==0) = 1;
E2 = E ./ colnorm;

% --- numeric rank (use dense SVD safely) ---
E2d   = full(E2);                    % svd does not support sparse
r_svd = rank(E2d);                   % numeric rank with default tol

% Singular values & condition estimate
[U,Sv,~] = svd(E2d, 'econ');
sv = diag(Sv);
sigma_max = sv(1);
sigma_min = sv(end);
cond_est  = sigma_max / max(sigma_min, eps);

fprintf('\n=== EFM basis diagnostics ===\n');
fprintf('Reactions (rows)    m = %d\n', m_eff);
fprintf('EFMs (cols, kept)   s = %d (dropped %d zero columns)\n', s_eff, sum(~col_nz));
fprintf('Numerical rank(E)   r = %d\n', r_svd);
fprintf('sigma_max(E)        = %.3g\n', sigma_max);
fprintf('sigma_min(E)        = %.3g\n', sigma_min);
fprintf('cond(E) ~           = %.3g\n', cond_est);
fprintf('Nullity (s - r)     = %d (redundant EFMs)\n', s_eff - r_svd);

tol_nz = 1e-12;                         % threshold for “active”
B = double(abs(E) > tol_nz);            % m x s binary (keeps sparsity)

% unique supports (exactly identical columns)
[~,~,ic] = unique(full(B).', 'rows', 'stable');   % unique needs full here
num_unique_supp = max(ic);

% structural rank of the pattern (not values)
if exist('sprank','file')
    r_struct = sprank(spones(B));       % best: maximum matching on sparsity pattern
else
    % fallback: QR on the 0/1 pattern
    [~,Rpat] = qr(spones(B), 0);
    r_struct = nnz(abs(diag(Rpat)) > 0);
end

fprintf('Unique supports     = %d / %d\n', num_unique_supp, s_eff);
fprintf('Rank(binary B)      = %d (structure-only)\n', r_struct);

% Per-anchor count 
if exist('EFM_anchor','var') && ~isempty(EFM_anchor)
    fprintf('Per-anchor accept counts (first 10 shown):\n');
    C = accumarray(EFM_anchor(:), 1);
    to_show = min(10, numel(C));
    disp(C(1:to_show).');
end

% Interpretation helper --> wil almost certainly always trigger a warning.
% EFMs like we enumerate them here are inherently collinear 
if r_svd <= 0.5*s_eff
    fprintf('Note: Many EFMs are numerically dependent; consider lowering JACCARD_TAU or\n');
    fprintf('      changing seeding to encourage more diverse supports.\n');
elseif sigma_min < 1e-10
    fprintf('Note: Basis is nearly singular; fitting may be ill-conditioned. Column scaling helps,\n');
    fprintf('      but you may need more/topologically distinct EFMs.\n');
end


%%
if exist('EFM_matrix','var') && ~isempty(EFM_matrix)
    tol_nz = eps_flux - tol_balance;          % same threshold you used elsewhere
    numEFMs = size(EFM_matrix, 2);

    % Active mask per reaction per EFM (exclude the residual row)
    active = abs(EFM_matrix(2:end, :)) > tol_nz;   % (n x numEFMs) logical

    % Build pair activity matrix: (p x numEFMs), true if either dir is active in that EFM
    p = numel(bpairs);
    pair_active = false(p, numEFMs);
    pair_names  = strings(p,1);
    for t = 1:p
        i = bpairs(t);
        j = badPartner(i);
        pair_active(t, :) = active(i, :) | active(j, :);
        pair_names(t) = sprintf('%s OR %s', rxnNames{i}, rxnNames{j});
    end

    % Pair-level coverage across all EFMs
    pair_covered = any(pair_active, 2);   % (p x 1)
    numCovered   = sum(pair_covered);
    fprintf('\nReversible-pair coverage: %d / %d pairs covered by >=1 EFM (%.1f%%)\n', ...
            numCovered, p, 100*numCovered/max(1,p));

    % List missing pairs (if any)
    if any(~pair_covered)
        missing_t = find(~pair_covered);
        fprintf('Pairs NOT covered by any EFM:\n');
        for t = reshape(missing_t,1,[])
            i = bpairs(t); j = badPartner(i);
            fprintf('  %s  OR  %s\n', rxnNames{i}, rxnNames{j});
        end
    else
        disp('All reversible pairs have at least one active direction in the EFM set.');
    end

   
end

%% 7c) Collapsed coverage: every original reaction covered by >=1 EFM?
% active: (n x numEFMs) logical, from earlier
numEFMs = size(EFM_matrix, 2);
tol_nz  = eps_flux - tol_balance;        % same threshold everywhere

% Per-reaction coverage across all EFMs
rxn_covered = any(active, 2);            % (n x 1)

% Reversible pairs
% bpairs: indices of the "first" member of each pair, p = numel(bpairs)
% badPartner(i): index of i's twin, or 0 if none
pair_names = strings(p,1);
for t = 1:p
    i = bpairs(t); j = badPartner(i);
    pair_names(t) = sprintf('%s OR %s', rxnNames{i}, rxnNames{j});
end

% Irreversibles (no twin)
irrevs = find(badPartner == 0);
numIr  = numel(irrevs);

% Build group-level activity matrix (groups = pairs first, then irrevs)
numGroups    = p + numIr;
group_names  = strings(numGroups, 1);
group_active = false(numGroups, numEFMs);

% Fill pairs
group_names(1:p)     = pair_names;
group_active(1:p, :) = pair_active;                  % (p x numEFMs), from earlier

% Fill irrevs
for s = 1:numIr
    idx = irrevs(s);
    g   = p + s;
    group_names(g)   = string(rxnNames{idx});
    group_active(g,:)= active(idx, :);               % (1 x numEFMs)
end

% Group coverage across the whole EFM set
group_covered = any(group_active, 2);                % (numGroups x 1)
numCovered    = sum(group_covered);

fprintf('\nCollapsed reaction coverage (pairs collapsed): %d / %d groups covered (%.1f%%)\n', ...
        numCovered, numGroups, 100*numCovered/max(1,numGroups));

if any(~group_covered)
    fprintf('Groups NOT covered by any EFM:\n');
    for g = find(~group_covered).'
        fprintf('  %s\n', group_names(g));
    end
else
    disp('All reactions are covered under the condition “either f or b active.”');
end


%%
idx = find(~covered);
fprintf('Uncovered reactions with GPR (%d):\n', numel(idx));
for t = 1:numel(idx)
    k = idx(t);
    rid = model_ir.rxns{k};
    gpr = '';
    if isfield(model_ir,'grRules') && numel(model_ir.grRules) >= k && ~isempty(model_ir.grRules{k})
        gpr = model_ir.grRules{k};
    elseif isfield(model_ir,'rxnGeneMat') && isfield(model_ir,'genes')
        genes = model_ir.genes(model_ir.rxnGeneMat(k,:) ~= 0);
        if ~isempty(genes), gpr = strjoin(genes,' OR '); end
    end
    if isempty(gpr), gpr = '(no GPR)'; end
    fprintf('%s | %s\n', rid, gpr);
end
%% — Helpers —
function x = getSolutionVector(sol)
    if isfield(sol,'x')   && ~isempty(sol.x)
        x = sol.x;
    elseif isfield(sol,'xn') && ~isempty(sol.xn)
        x = sol.xn;
    elseif isfield(sol,'pool') && ~isempty(sol.pool) && isfield(sol.pool(1),'xn')
        x = sol.pool(1).xn;
    else
        error('No solution vector found.');
    end
end

function T = applyNoGoodCut(T, supp, n)
    oth = setdiff(1:n, supp);
    cut = sparse(1,2*n);
    cut(n+supp) = 1;
    cut(n+oth)  = -1;
    T.A     = [T.A; cut];
    T.rhs   = [T.rhs; numel(supp)-1];
    T.sense = [T.sense; '<'];
end

function tf = no_bad_pairs(supp, badPartner)
% true if no {i,j} with j = badPartner(i) inside supp
    tf = true;
    if isempty(supp), return; end
    bp = badPartner(supp);
    if any(bp>0)
        tf = ~any(ismember(bp(bp>0), supp));
    end
end

function supp_out = simplify_bad_pairs(supp_in, badPartner, keepIdx)
% Remove one from each conflicting pair, prefer to keep 'keepIdx' if involved
    supp = supp_in(:)';
    changed = true;
    while changed
        changed = false;
        for t = 1:numel(supp)
            i = supp(t);
            j = badPartner(i);
            if j>0 && ismember(j, supp)
                % both i and j in supp → drop one
                if i == keepIdx
                    drop = j;
                elseif j == keepIdx
                    drop = i;
                else
                    % arbitrary but stable: drop the larger index
                    drop = max(i,j);
                end
                supp(supp==drop) = [];
                changed = true;
                break;
            end
        end
    end
    supp_out = supp(:);
end

function ok = nullity1_ok(S, supp)
    if isempty(supp), ok=false; return; end
    [~,R] = qr(S(:,supp),0);         % sparse QR
    rk = sum(abs(diag(R)) > 1e-10);  % robust tol
    ok = (numel(supp) - rk == 1);
end
function [ok, supp_out, v_full] = prune_nullity_by_single_drop( ...
            S, ub_eff, supp0, eps_flux, tol_balance, badPartner, anchor)

    ok = false; supp_out = []; v_full = [];

    if numel(supp0) <= 1, return; end

    
    order = supp0(:).';
    if nargin >= 7 && ~isempty(anchor)
        order = [ order(order ~= anchor), anchor ];
    end

    for j = order
        if nargin >= 7 && ~isempty(anchor) && j == anchor
            % Prefer to keep anchor; only drop if nothing else works
            if numel(supp0) > 2
                % we'll try anchor last anyway
            end
        end

        trial = setdiff(supp0, j, 'stable');
        if isempty(trial), continue; end

        % no bad pairs on trial
        if ~no_bad_pairs(trial, badPartner), continue; end

        % anchor must remain
        if nargin >= 7 && ~isempty(anchor) && ~ismember(anchor, trial)
            continue;
        end

        % LP on trial: Sv=0, eps<=v<=U
        [ok_refit, v_trial] = refit_on_support_feasible(S, ub_eff, trial, eps_flux);
        if ~ok_refit, continue; end

        % conservative nullity on trial
        if ~nullity1_ok(S, trial), continue; end

        % assemble full vector & check FULL Sv
        n = size(S,2);
        v_full = zeros(n,1);
        v_full(trial) = v_trial;
        if norm(S(:,trial)*v_full(trial), Inf) > 1e-5 %tol_balance
            disp("dang")
            continue;
        end

        ok = true; supp_out = trial(:); return;
    end
end

function [ok, v_s] = refit_on_support_feasible(S, ub_eff, supp, eps_pos)
% Solve: S(:,supp)*v = 0,  eps_pos <= v <= Uj  (no sum constraint)
    m = size(S,1);
    k = numel(supp);
    if k == 0, ok=false; v_s=[]; return; end

    model.A      = sparse(S(:,supp));
    model.rhs    = zeros(m,1);
    model.sense  = repmat('=', m, 1);
    model.lb     = eps_pos * ones(k,1);
    model.ub     = ub_eff(supp);
    model.vtype  = repmat('C', k, 1);
    model.modelsense = 'min';
    model.obj    = zeros(k,1);         % any feasible solution

    sol = gurobi(model, struct('OutputFlag',0,'FeasibilityTol',1e-9));
    ok  = isfield(sol,'status') && strcmp(sol.status,'OPTIMAL');
    if ok, v_s = sol.x; else, v_s = []; end
end

function [res_supp, maxOff, idxOff] = support_debug_metrics(S, v_raw, supp)
    if isempty(supp)
        res_supp = 0;
    else
        res_supp = norm(S(:,supp) * v_raw(supp), Inf);
    end
    n = size(S,2);
    off = setdiff(1:n, supp);

    if isempty(off)
        maxOff = 0; idxOff = [];
        return;
    end

    voff = abs(v_raw(off));
    [maxOff, k] = max(voff);
    idxOff = off(k);
    
end


function [supp_out, v_new] = cancel_twins(supp_in, v, badPartner, S, eps, tol)
    v_new = v(:);
    supp  = supp_in(:).';
    n = numel(v_new);
    keepWinner = false(n,1);
    seen = false(n,1);

    for t = 1:numel(supp)
        i = supp(t); j = badPartner(i);
        if j<=0 || ~ismember(j, supp) || seen(i) || seen(j), continue; end
        seen([i j]) = true;

        % only cancel if columns are (numerically) opposite
        if max(abs(S(:,i) + S(:,j))) <= 1e-9
            vi = v_new(i); vj = v_new(j);
            if vi > 0 && vj > 0
                if vi == vj
                    keepWinner(i) = true;
                    keepWinner(j) = true;
                    %disp(vi)
                    %disp(vj)
                    continue;
                end 

                d = min(vi, vj);
                v_new(i) = vi - d;
                v_new(j) = vj - d;

                % winner has the residual > 0; keep it in support even if < eps
                if v_new(i) > v_new(j) && v_new(i) > tol
                    keepWinner(i) = true;
                elseif v_new(j) > tol
                    keepWinner(j) = true;
                end
            end
        end
    end

   
    supp_out = unique([ find(v_new >= eps - tol); find(keepWinner) ]);
end


function [hasAny, idxPairs] = badpairs_in_support(supp, badPartner, bpairs)
% Return whether the support has any twin pair, and which bpairs indices
    if isempty(supp)
        hasAny = false; idxPairs = []; return;
    end
    inSupp = false(size(badPartner));
    inSupp(supp) = true;
    j = badPartner(bpairs);
    present = inSupp(bpairs) & inSupp(j);
    hasAny = any(present);
    idxPairs = find(present);
end


function [ok, v_full, objval] = maximize_sum_v_on_support(S, ub_eff, supp, eps_pos)
% Maximize total flux on a fixed support:
%   maximize sum(v_supp)  s.t.  S(:,supp)*v = 0,  eps_pos <= v <= ub_eff(supp)
    m = size(S,1); n = size(S,2);
    k = numel(supp);
    if k==0, ok=false; v_full=[]; objval=-Inf; return; end

    model.A          = sparse(S(:,supp));
    model.rhs        = zeros(m,1);
    model.sense      = repmat('=', m, 1);
    model.lb         = eps_pos*ones(k,1);
    model.ub         = ub_eff(supp);
    model.vtype      = repmat('C', k, 1);
    model.modelsense = 'max';
    model.obj        = ones(k,1);           % maximize sum of fluxes

    params = struct('OutputFlag',0,'FeasibilityTol',1e-9);
    sol = gurobi(model, params);
    ok  = isfield(sol,'status') && strcmp(sol.status,'OPTIMAL');
    if ok
        v_full = zeros(n,1);
        v_full(supp) = sol.x;
        objval = sol.objval;
    else
        v_full = []; objval = -Inf;
    end
end


function J = jaccard_max_sim(supp, accepted_supps)
% Max Jaccard similarity between 'supp' and any support in 'accepted_supps'
    J = 0;
    A = unique(supp(:));
    for t = 1:numel(accepted_supps)
        B = unique(accepted_supps{t}(:));
        inter = numel(intersect(A,B));
        uni   = numel(union(A,B));
        j = inter / max(uni,1);
        if j > J, J = j; end
        if J >= 0.999, return; end
    end
end


function save_checkpoint_atomic(targetFile, model_sig)
% Write state to a temp file, then atomically rename → targetFile.
    try rng_state = rng; catch, rng_state = []; end 
    model_sig = model_sig; 

    % Variables to persist 
    varNames = { ...
        'EFM_supps','EFM_vecs','EFM_anchor','accepted_supps', ...
        'covered','skipped','linModel', ...
        'noProgressCount','prevCoveredCount','stallTerminated', ...
        'efm_with_badpair_cnt','badpair_counts', ...
        'JACCARD_TAU','MAX_PER_ANCHOR','MIN_PER_ANCHOR', ...
        'eps_flux','tol_balance','M', ...
        'ub','badPartner','bpairs','p', ...
        'S','rxnNames','model_sig','rng_state'};

    Sstruct = struct();
    for v = varNames
        vn = v{1};
        if evalin('caller', sprintf('exist(''%s'',''var'')', vn))
            Sstruct.(vn) = evalin('caller', vn);
        end
    end

    tmp = [targetFile, '.tmp'];
    save(tmp, '-struct', 'Sstruct', '-v7.3');   
    movefile(tmp, targetFile, 'f');             
    fprintf('[checkpoint] %s\n', targetFile)
en

function load_checkpoint_state(fname)
% Load all fields from checkpoint into caller workspace.
    L = load(fname);
    f = fieldnames(L);
    for k = 1:numel(f), assignin('caller', f{k}, L.(f{k})); end
    fprintf('[resume] state loaded from %s\n', fname);
end
