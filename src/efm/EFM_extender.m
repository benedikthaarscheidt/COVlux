
% EFM extender with proper EFM loading and random anchor selection
% Based on proven EFM loading method from working script

%% 0) Gurobi threads & paths
nThreads = feature('numcores')-2;
fprintf('Using %d threads for Gurobi\n',nThreads);
setenv('OMP_NUM_THREADS',num2str(nThreads));
addpath(fullfile(getenv('GUROBI_HOME'),'examples','matlab'),'-begin');

%% 1) Load irreversible pruned model
data     = load('/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/iML1515_pruned_permissive_biomassrm_loopless_v5.mat','pruned_ir');
model_ir = data.pruned_ir;
S        = model_ir.S;
rxnNames = model_ir.rxns;
[m,n]    = size(S);

%% 2) Load existing EFMs PROPERLY (based on working script)
efmMatPath1 = '/Users/benedikthaarscheidt/M.Sc./master_thesis/efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5_higher_cover.mat';
efmMatPath2 = '/Users/benedikthaarscheidt/M.Sc./master_thesis/scripts/EFM/efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5.mat';

fprintf('Loading EFMs using proven method...\n');

% Load both files exactly as in working script
S1 = load(efmMatPath1, 'EFM_matrix', 'rxnNames');
S2 = load(efmMatPath2, 'EFM_matrix', 'rxnNames');

% Extract EFM matrices and reaction names
E1 = S1.EFM_matrix(2:end, :);  % Remove residual row
rxn1 = string(S1.rxnNames(:));
E2 = S2.EFM_matrix(2:end, :);  
rxn2 = string(S2.rxnNames(:));

% Find intersection of reactions and combine EFMs
[rxnE, ia1, ia2] = intersect(rxn1, rxn2, 'stable');
E_full = [E1(ia1, :), E2(ia2, :)];

[mE, s_full] = size(E_full);
fprintf('Combined EFM matrix: %d reactions x %d EFMs\n', mE, s_full);

%% Extract supports from the combined EFM matrix
tol_nz = 1e-6; % Standard threshold
all_known_supps = cell(1, s_full);
all_known_vecs = cell(1, s_full);

for c = 1:s_full
    v = E_full(:, c);
    supp = find(abs(v) > tol_nz);
    all_known_supps{c} = supp(:);
    all_known_vecs{c} = v(:);
end

fprintf('Total loaded EFM supports: %d\n', numel(all_known_supps));

% Remove duplicate supports
unique_supps = {};
unique_vecs = {};
supp_hashes = {};

for k = 1:numel(all_known_supps)
    supp_k = sort(all_known_supps{k}(:));
    supp_hash = sprintf('%d,', supp_k);
    
    if ~ismember(supp_hash, supp_hashes)
        unique_supps{end+1} = supp_k;
        unique_vecs{end+1} = all_known_vecs{k};
        supp_hashes{end+1} = supp_hash;
    end
end

fprintf('After deduplication: %d unique EFM supports\n', numel(unique_supps));

% Calculate coverage
existing_coverage = false(mE, 1);
for k = 1:numel(unique_supps)
    existing_coverage(unique_supps{k}) = true;
end
fprintf('Coverage from existing EFMs: %d/%d reactions covered\n', sum(existing_coverage), mE);

%% 3) Parameters
eps_flux      = 1e-5;
tol_balance   = 1e-6;
M             = 1e3;
badPair_count = 0;
JACCARD_TAU   = 0.70;              % Lower threshold for more diversity
MAX_PER_ANCHOR = 150;              % More per anchor for new EFMs
MIN_PER_ANCHOR = 2;

accepted_supps = unique_supps;     % Start with known supports
EFM_supps  = unique_supps;         % all accepted supports (any anchor)
EFM_vecs   = unique_vecs;          % all accepted vectors
EFM_anchor = -ones(size(unique_supps)); % -1 for existing EFMs

%% Compute range of all reactions for fixing the flux 
fprintf('Computing reaction flux ranges...\n');
fvaModel.A        = sparse(S);
fvaModel.rhs      = zeros(m,1);
fvaModel.sense    = repmat('=', m, 1);
fvaModel.lb       = model_ir.lb;
fvaModel.ub       = model_ir.ub;
fvaModel.vtype    = repmat('C', n, 1);
fvaModel.modelsense = 'min';
fvaModel.obj      = zeros(n,1);

% Gurobi parameters for FVA
fvaParams.OutputFlag = 0;
fvaParams.TimeLimit  = 60;
fvaParams.FeasibilityTol = 1e-9;

vMin = nan(n,1);
vMax = nan(n,1);

for j = 1:n
    % minimize v_j
    fvaModel.obj(:) = 0;
    fvaModel.obj(j) = 1;
    fvaModel.modelsense = 'min';
    sol = gurobi(fvaModel, fvaParams);
    if isfield(sol,'status') && strcmp(sol.status,'OPTIMAL')
        vMin(j) = sol.objval;
    else
        warning('FVA min infeasible for reaction %d (%s)', j, rxnNames{j});
    end

    % maximize v_j
    fvaModel.modelsense = 'max';
    sol = gurobi(fvaModel, fvaParams);
    if isfield(sol,'status') && strcmp(sol.status,'OPTIMAL')
        vMax(j) = sol.objval;
    else
        warning('FVA max infeasible for reaction %d (%s)', j, rxnNames{j});
    end
end

range = [vMin, vMax];

% Handle infinite bounds
vMin(~isfinite(vMin)) = 0;
vMax(~isfinite(vMax)) = M;
vMin = max(vMin, 0);
range = [vMin, vMax];

%% 4) Build static MILP template 
fprintf('Building MILP template...\n');

lb = [range(:,1); zeros(n,1)];
ub = [range(:,2);  ones(n,1)];
ub_eff = ub; 
ub_eff(~isfinite(ub_eff)) = M;
ub = ub_eff;

vtype = [repmat('C',n,1); repmat('B',n,1)];

Aeq = [sparse(S), sparse(m,n)];
beq = zeros(m,1);

% v_i - M*b_i <= 0
% -v_i + eps*b_i <= 0
A1 = [speye(n), -ub(1:n).*speye(n);
     -speye(n), eps_flux*speye(n)];
b1 = zeros(2*n,1);

% Import/export constraints
nz = (S~=0);
importCols = false(1,n);
exportCols = false(1,n);
for j = 1:n
    rows = nz(:,j);
    if any(rows)
        col = S(rows,j);
        importCols(j) = all(col>0);
        exportCols(j) = all(col<0);
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

% Bad pair constraints
badPartner = zeros(n,1);
badPair_count = 0;
for i = 1:n
    rn = rxnNames{i};
    if endsWith(rn,'_f')
        j = find(strcmp(rxnNames,[rn(1:end-2) '_b']),1);
        if ~isempty(j)
            badPartner(i)=j;
            badPartner(j)=i;
            if i < j
                badPair_count = badPair_count + 1;
            end
        end
    end
end

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

% Combine all constraints
Aineq = [A1; Aimp; Aexp; Abp];
bineq = [b1; bimp; bexp; bbp];

sense = [repmat('=', size(Aeq,1),1); repmat('<', size(Aineq,1),1)];

obj = [zeros(n,1); ones(n,1)];

% Build base MILP model
linModel.A          = [Aeq; Aineq];
linModel.rhs        = [beq; bineq];
linModel.sense      = sense;
linModel.lb         = lb;
linModel.ub         = ub;
linModel.vtype      = vtype;
linModel.modelsense = 'min';
linModel.obj        = obj;

% Add no-good cuts for ALL known supports
fprintf('Adding no-good cuts for %d known EFM supports...\n', numel(unique_supps));
for k = 1:numel(unique_supps)
    linModel = applyNoGoodCut(linModel, unique_supps{k}, n);
end
fprintf('Added %d no-good cuts\n', numel(unique_supps));

%% 5) Gurobi parameters
commonParams = struct( ...
   'FeasibilityTol',1e-9, ...
   'IntFeasTol',    1e-9, ...
   'NumericFocus',  3);

pSeed = struct( ...
  'OutputFlag',     0, ...
  'PoolSearchMode', 2, ...
  'PoolSolutions',  100, ...  % More candidates
  'TimeLimit',      300, ...  % Longer time
  'MIPGap',         0.3, ...  % Accept suboptimal
  'Cuts',           -1, ... 
  'Presolve',       2, ...
  'Heuristics',     0.8, ...  % More aggressive
  'MIPFocus',       1, ...
  'Threads',        nThreads);

for p = {pSeed}
    p{1}.FeasibilityTol = commonParams.FeasibilityTol;
    p{1}.IntFeasTol     = commonParams.IntFeasTol;
    p{1}.NumericFocus   = commonParams.NumericFocus;
end

%% 6) Checkpoint setup
CHECKPOINT_ENABLE  = true;
AUTOSAVE_SEC       = 15*60;  % Save every 15 minutes
CHECKPOINT_FILE    = fullfile(pwd,'efm_checkpoint_NEW.mat');

% Model signature for resume validation
model_sig = struct('m',m,'n',n,'nnzS',nnz(S),'rxnCount',numel(rxnNames));

% Optional: resume from checkpoint
LOAD_CHECKPOINT = true;
LOAD_FILE       = CHECKPOINT_FILE;

if LOAD_CHECKPOINT && exist(LOAD_FILE,'file')
    R = load(LOAD_FILE,'model_sig');
    if isfield(R,'model_sig') && isequal(R.model_sig, model_sig)
        fprintf('[resume] Loading: %s\n', LOAD_FILE);
        load_checkpoint_state(LOAD_FILE); 
        if exist('rng_state','var'), rng(rng_state); end
    else
        warning('[resume] Incompatible checkpoint; starting fresh.');
    end
end

% Initialize tracking variables if not resuming
if ~exist('covered','var'), covered = existing_coverage; end
if ~exist('skipped','var'), skipped = false(mE, 1); end
if ~exist('noProgressCount','var'), noProgressCount = 0; end
if ~exist('prevCoveredCount','var'), prevCoveredCount = sum(covered); end

lastSave = tic;

%% 7) RANDOM anchor selection enumeration
fprintf('\n=== STARTING RANDOM ANCHOR SEARCH ===\n');

stallLimit = 500;
stallTerminated = false;
maxseedattempts=20;

% Search all reactions for novel combinations
all_searchable_rxns = 1:mE;
fprintf('Searching across all %d reactions for novel EFM patterns\n', mE);

try
    while any(~skipped)
        % Checkpoint save
        if CHECKPOINT_ENABLE && toc(lastSave) > AUTOSAVE_SEC
            save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
            lastSave = tic;
        end

        % Random selection from all searchable reactions
        active_rxns = all_searchable_rxns(~skipped);
        if isempty(active_rxns)
            break;
        end
        
        idx = randi(numel(active_rxns));
        i = active_rxns(idx);
        fprintf('\nTrying random anchor: %d (%s) - Coverage: %d/%d\n', ...
                i, rxnE{i}, sum(covered), mE);
        
        localTemplate = linModel;
        found = false;
        accepted_this_anchor = 0;

        for attempt = 1:maxSeedAttempts
            M1 = localTemplate;
            M1.lb(n+i) = 1;
            M1.ub(n+i) = 1;

            sol1 = gurobi(M1, pSeed);
            if ~isfield(sol1,'pool') || isempty(sol1.pool)
                fprintf('  attempt %d: no pool solutions → retrying\n', attempt);
                continue;
            end

            numCand = numel(sol1.pool);
            fprintf('  attempt %d: %d pool candidates\n', attempt, numCand);

            supp_list = cell(numCand,1);
            v_list = cell(numCand,1);
            hasPair_list = false(numCand,1);
            pairIdx_list = cell(numCand,1);
            badpair_violate = false(numCand,1);
            acceptType = zeros(numCand,1);
            supp_accept = cell(numCand,1);
            v_accept = cell(numCand,1);
            anchorInSupp = false(numCand,1);

            parfor k = 1:numCand
                xk = getSolutionVector(sol1.pool(k));
                v_k = xk(1:n);
                supp0 = find(v_k >= eps_flux - tol_balance);
                if isempty(supp0) || ~ismember(i, supp0)
                    continue;
                end
                if ~no_bad_pairs(supp0, badPartner)
                    [supp0, v_k] = cancel_twins(supp0, v_k, badPartner, S, eps_flux, tol_balance);
                    if isempty(supp0) || ~ismember(i, supp0)
                        continue;
                    end
                end
                anchorInSupp(k) = true;
                [hBP, pIdx] = badpairs_in_support(supp0, badPartner, bpairs);
                hasPair_list(k) = hBP;
                if hBP, pairIdx_list{k} = pIdx; end

                nullOK = nullity1_ok(S, supp0);
                res_supp = norm(S(:,supp0) * v_k(supp0), Inf);
                full_suppOK = (res_supp <= 1e-5);
                stillBP = ~no_bad_pairs(supp0, badPartner);
                badpair_violate(k) = stillBP;

                if nullOK && full_suppOK && ~stillBP
                    v_full = v_k;
                    acceptType(k) = 1;
                    supp_accept{k} = supp0;
                    v_accept{k} = v_full;
                    continue;
                end

                if (~nullOK || ~full_suppOK)
                    [ok_fix, supp_fix, v_fix] = prune_nullity_by_single_drop(S, ub, supp0, eps_flux, tol_balance, badPartner, i);
                    if ok_fix
                        res_supp_fix = norm(S(:,supp_fix) * v_fix(supp_fix), Inf);
                        full_suppOK_fix = (res_supp_fix <= 1e-5);
                        stillBP_fix = ~no_bad_pairs(supp_fix, badPartner);
                        if full_suppOK_fix && ~stillBP_fix
                            acceptType(k) = 2;
                            supp_accept{k} = supp_fix(:);
                            v_accept{k} = v_fix(:);
                        end
                    end
                end
            end

            for k = 1:numCand
                if ~anchorInSupp(k), continue; end

                if badpair_violate(k)
                    continue;
                end

                if acceptType(k) == 1 || acceptType(k) == 2
                    if acceptType(k) == 1
                        supp0 = supp_accept{k};
                        v_k = v_accept{k};
                        tag = '';
                    else
                        supp0 = supp_accept{k};
                        v_k = v_accept{k};
                        tag = ' (single-drop)';
                    end

                    % Jaccard similarity screen
                    J = 0;
                    if ~isempty(accepted_supps)
                        J = jaccard_max_sim(supp0, accepted_supps);
                    end
                    if J >= JACCARD_TAU
                        localTemplate = applyNoGoodCut(localTemplate, supp0, n);
                        continue;
                    end

                    % Use Gurobi solution directly (no expensive refit)
                    v_full = zeros(n,1);
                    v_full(supp0) = v_k(supp0);

                    % Record globally
                    EFM_supps{end+1} = supp0(:);
                    EFM_vecs{end+1} = v_full(:);
                    EFM_anchor(end+1) = i;
                    accepted_supps{end+1} = supp0(:);

                    % Update coverage and forbid support
                    covered(supp0) = true;
                    localTemplate = applyNoGoodCut(localTemplate, supp0, n);

                    fprintf('   → NEW EFM%s |supp|=%d | Sv=%.2e | sum v=%.3g | J=%.2f\n', ...
                            tag, numel(supp0), norm(S*v_full,Inf), sum(v_full(supp0)), J);

                    accepted_this_anchor = accepted_this_anchor + 1;

                    % Checkpoint after acceptance
                    if CHECKPOINT_ENABLE && toc(lastSave) > AUTOSAVE_SEC
                        save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
                        lastSave = tic;
                    end
                else
                    % Reject - add no-good cut
                    localTemplate = applyNoGoodCut(localTemplate, supp_list{k}, n);
                end
            end

            % Stop if we hit the cap
            if accepted_this_anchor >= MIN_PER_ANCHOR
                break;
            end
        end

        % Mark as skipped if nothing found
        if accepted_this_anchor == 0
            fprintf('   No valid NEW support for reaction %d after %d attempts\n', i, maxSeedAttempts);
            skipped(i) = true;
        end

        % Progress tracking
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
            
            % Final checkpoint before breaking
            if CHECKPOINT_ENABLE
                save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
            end
            break;
        end
    end
catch ME
    fprintf('[error] %s\n', ME.message);
    if CHECKPOINT_ENABLE
        try
            save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
        catch
        end
    end
    rethrow(ME);
end

% Final checkpoint
if CHECKPOINT_ENABLE
    save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
end

%% 8) Build complete results
fprintf('\nBuilding complete EFM matrix...\n');

% Convert to proper EFM matrix format
EFM_matrix = zeros(mE+1, numel(EFM_supps));
EFM_support = false(mE, numel(EFM_supps));

for c = 1:numel(EFM_supps)
    v = EFM_vecs{c};
    EFM_matrix(2:end, c) = v;
    EFM_matrix(1, c) = norm(S*v, Inf);
    EFM_support(:, c) = (abs(v) > tol_nz);
end

% Create names
rowNames = [{'Sv_inf_residual'}; rxnE(:)];
varNames = cell(1, numel(EFM_supps));
for c = 1:numel(EFM_supps)
    if EFM_anchor(c) == -1
        varNames{c} = sprintf('EXISTING_EFM_%d', c);
    else
        varNames{c} = sprintf('EFM_%d_anchor_%s', c, rxnE{EFM_anchor(c)});
    end
end
varNames = matlab.lang.makeValidName(varNames, 'ReplacementStyle','delete');
varNames = matlab.lang.makeUniqueStrings(varNames);

EFM_table = array2table(EFM_matrix, 'RowNames', rowNames, 'VariableNames', varNames);

%% 9) Report & save results
numExisting = numel(unique_supps);
numNew = numel(EFM_supps) - numExisting;

fprintf('\n=== FINAL RESULTS ===\n');
fprintf('Existing EFMs: %d\n', numExisting);
fprintf('New EFMs found: %d\n', numNew);
fprintf('Total EFMs: %d\n', numel(EFM_supps));
fprintf('Final coverage: %d/%d reactions\n', sum(covered), mE);

% Save combined set
save('efms_COMBINED_iML1515_pruned_permissive_biomassrm_loopless_v5.mat', ...
     'EFM_matrix', 'EFM_table', 'rowNames', 'varNames', 'EFM_support', ...
     'EFM_anchor', 'EFM_supps', 'EFM_vecs', 'rxnE', 'S', 'model_ir', '-v7.3');

fprintf('Saved combined set (%d EFMs)\n', numel(EFM_supps));

%% 

% Save new EFMs separately if any found
if numNew > 0
    new_idx = numExisting+1:numel(EFM_supps);
    new_EFM_matrix = EFM_matrix(:, new_idx);
    new_EFM_table = EFM_table(:, new_idx);
    new_EFM_support = EFM_support(:, new_idx);
    new_EFM_anchor = EFM_anchor(new_idx);
    new_EFM_supps = EFM_supps(new_idx);
    new_EFM_vecs = EFM_vecs(new_idx);
    new_varNames = varNames(new_idx);

    new_coverage = false(mE, 1);
    for k = 1:numel(new_EFM_supps) % Loop through just the new supports
        supp = new_EFM_supps{k};
        if ~isempty(supp)
            new_coverage(supp) = true;
        end
    end
    new_coverage_count = sum(new_coverage);
    
    fprintf('Coverage from NEW EFMs only: %d/%d reactions\n', new_coverage_count, mE);

    
    save('efms_NEW_iML1515_pruned_permissive_biomassrm_loopless_v5.mat', ...
         'new_EFM_matrix', 'new_EFM_table', 'rowNames', 'new_varNames', ...
         'new_EFM_support', 'new_EFM_anchor', 'new_EFM_supps', 'new_EFM_vecs', ...
         'new_coverage', 'new_coverage_count', ... % <-- Added new variables
         'rxnE', 'S', 'model_ir', '-v7.3');
    
    fprintf('Saved %d NEW EFMs to separate file\n', numNew);
end
fprintf('\n=== COMPLETE ===\n');

%% Helper functions
function x = getSolutionVector(sol)
    if isfield(sol,'x') && ~isempty(sol.x)
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
    cut(n+oth) = -1;
    T.A = [T.A; cut];
    T.rhs = [T.rhs; numel(supp)-1];
    T.sense = [T.sense; '<'];
end

function tf = no_bad_pairs(supp, badPartner)
    tf = true;
    if isempty(supp), return; end
    bp = badPartner(supp);
    if any(bp>0)
        tf = ~any(ismember(bp(bp>0), supp));
    end
end

function ok = nullity1_ok(S, supp)
    if isempty(supp), ok=false; return; end
    [~,R] = qr(S(:,supp),0);
    rk = sum(abs(diag(R)) > 1e-10);
    ok = (numel(supp) - rk == 1);
end

function [ok, supp_out, v_full] = prune_nullity_by_single_drop(S, ub_eff, supp0, eps_flux, tol_balance, badPartner, anchor)
    ok = false; supp_out = []; v_full = [];
    if numel(supp0) <= 1, return; end

    order = supp0(:).';
    if nargin >= 7 && ~isempty(anchor)
        order = [order(order ~= anchor), anchor];
    end

    for j = order
        if nargin >= 7 && ~isempty(anchor) && j == anchor
            if numel(supp0) > 2
                continue;
            end
        end

        trial = setdiff(supp0, j, 'stable');
        if isempty(trial), continue; end
        if ~no_bad_pairs(trial, badPartner), continue; end
        if nargin >= 7 && ~isempty(anchor) && ~ismember(anchor, trial)
            continue;
        end

        [ok_refit, v_trial] = refit_on_support_feasible(S, ub_eff, trial, eps_flux);
        if ~ok_refit, continue; end
        if ~nullity1_ok(S, trial), continue; end

        n = size(S,2);
        v_full = zeros(n,1);
        v_full(trial) = v_trial;
        if norm(S(:,trial)*v_full(trial), Inf) > 1e-5
            continue;
        end

        ok = true; supp_out = trial(:); return;
    end
end

function [ok, v_s] = refit_on_support_feasible(S, ub_eff, supp, eps_pos)
    m = size(S,1);
    k = numel(supp);
    if k == 0, ok=false; v_s=[]; return; end

    model.A = sparse(S(:,supp));
    model.rhs = zeros(m,1);
    model.sense = repmat('=', m, 1);
    model.lb = eps_pos * ones(k,1);
    model.ub = ub_eff(supp);
    model.vtype = repmat('C', k, 1);
    model.modelsense = 'min';
    model.obj = zeros(k,1);

    sol = gurobi(model, struct('OutputFlag',0,'FeasibilityTol',1e-9));
    ok = isfield(sol,'status') && strcmp(sol.status,'OPTIMAL');
    if ok, v_s = sol.x; else, v_s = []; end
end

function [supp_out, v_new] = cancel_twins(supp_in, v, badPartner, S, eps, tol)
    v_new = v(:);
    supp = supp_in(:).';
    n = numel(v_new);
    keepWinner = false(n,1);
    seen = false(n,1);

    for t = 1:numel(supp)
        i = supp(t); j = badPartner(i);
        if j<=0 || ~ismember(j, supp) || seen(i) || seen(j), continue; end
        seen([i j]) = true;

        if max(abs(S(:,i) + S(:,j))) <= 1e-9
            vi = v_new(i); vj = v_new(j);
            if vi > 0 && vj > 0
                if vi == vj
                    keepWinner(i) = true;
                    keepWinner(j) = true;
                    continue;
                end

                d = min(vi, vj);
                v_new(i) = vi - d;
                v_new(j) = vj - d;

                if v_new(i) > v_new(j) && v_new(i) > tol
                    keepWinner(i) = true;
                elseif v_new(j) > tol
                    keepWinner(j) = true;
                end
            end
        end
    end

    supp_out = unique([find(v_new >= eps - tol); find(keepWinner)]);
end

function [hasAny, idxPairs] = badpairs_in_support(supp, badPartner, bpairs)
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

function J = jaccard_max_sim(supp, accepted_supps)
    J = 0;
    A = unique(supp(:));
    for t = 1:numel(accepted_supps)
        B = unique(accepted_supps{t}(:));
        inter = numel(intersect(A,B));
        uni = numel(union(A,B));
        j = inter / max(uni,1);
        if j > J, J = j; end
        if J >= 0.999, return; end
    end
end

%% Checkpoint functions
function save_checkpoint_atomic(targetFile, model_sig)
    try 
        rng_state = rng; 
    catch
        rng_state = []; 
    end

    varNames = { ...
        'EFM_supps','EFM_vecs','EFM_anchor','accepted_supps', ...
        'covered','skipped','linModel', ...
        'noProgressCount','prevCoveredCount','stallTerminated', ...
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
    fprintf('[checkpoint] %s\n', targetFile);
end

function load_checkpoint_state(fname)
    L = load(fname);
    f = fieldnames(L);
    for k = 1:numel(f)
        assignin('caller', f{k}, L.(f{k}));
    end
    fprintf('[resume] state loaded from %s\n', fname);
end