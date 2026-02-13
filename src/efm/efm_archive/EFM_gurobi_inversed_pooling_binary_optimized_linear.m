% with immediate no‐good cuts for invalid supports only,
% vector bad‐pair test, seeding with retries, pooling without cardinality cuts

%% 0) Gurobi threads & paths
nThreads = feature('numcores');
fprintf('Using %d threads for Gurobi\n',nThreads);
setenv('OMP_NUM_THREADS',num2str(nThreads));
addpath(fullfile(getenv('GUROBI_HOME'),'examples','matlab'),'-begin');

%% 1) Load irreversible pruned model
data     = load('/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/iML1515_pruned_permissive_biomassrm_loopless_v5.mat','pruned_ir')
model_ir = data.pruned_ir;
S        = model_ir.S;
rxnNames = model_ir.rxns;
[m,n]    = size(S);

%% 2) Parameters
eps_flux      = 1e-5;
tol_balance   = 1e-6;
M             = 1e3;
badPair_count = 0;

%% Compute range of all reactions for fixing the flux 

fvaModel.A        = sparse(S);
fvaModel.rhs      = zeros(m,1);
fvaModel.sense    = repmat('=', m, 1);
fvaModel.lb       = model_ir.lb;
fvaModel.ub       = model_ir.ub;
fvaModel.vtype    = repmat('C', n, 1);
fvaModel.modelsense = 'min';      % we'll switch between 'min' and 'max'
fvaModel.obj      = zeros(n,1);

% Gurobi parameters for FVA
fvaParams.OutputFlag = 0;
fvaParams.TimeLimit  = 60;    % e.g. 60 s per LP
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

% Optional: show a table with vMin/vMax for these
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
pSeed = struct( ...
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

maxSeedAttempts = 20;
stallLimit = 500;
noProgressCount   = 0;
prevCoveredCount  = sum(covered);
stallTerminated   = false;


pool = gcp('nocreate');
if isempty(pool), pool = parpool('local', nThreads); end

while any(~covered & ~skipped)
    toTry = find(~covered & ~skipped);
    fprintf('\n=== Next reaction: %d remaining (%d covered) ===\n', numel(toTry), sum(covered));

    % pick next uncovered reaction
    idx = randi(numel(toTry));
    i = toTry(idx);
    prevUc = i; 

    localTemplate = linModel;
    found = false;

    for attempt = 1:maxSeedAttempts
        % force b_i = 1
        M1 = localTemplate;
        M1.lb(n+i) = 1;  
        M1.ub(n+i) = 1;

        sol1 = gurobi(M1, pSeed);
        if ~isfield(sol1,'pool') || isempty(sol1.pool)
            fprintf('  attempt %d: no pool solutions → retrying\n', attempt);
            continue;
        end

        numCand  = numel(sol1.pool);
        fprintf('  attempt %d: %d pool candidates\n', attempt, numCand);

        accepted_this_attempt = false;

        for k = 1:numCand
            xk    = getSolutionVector(sol1.pool(k));
            v_k   = xk(1:n);
            b   = xk(n+1:end);
            supp0 = find(v_k >= eps_flux-tol_balance);   
            %supp0 = find(b>0.5);   

            if isempty(supp0) || ~ismember(i, supp0), continue; end
            
            % Quick bad-pair screen (prefer keeping i)
            if ~no_bad_pairs(supp0, badPartner)
                [supp0, v_k] = cancel_twins(supp0, v_k, badPartner, S, eps_flux, tol_balance);
            end
            
            % Accept-if-fits: strict nullity + FULL steady state
            nullOK   = nullity1_ok(S, supp0);
            res_full = norm(S * v_k, Inf);
            fullOK   = (res_full <= tol_balance);
            res_supp= norm(S(:,supp0) * v_k(supp0), Inf);
            full_suppOK=(res_supp <= 1e-5);
            if nullOK && full_suppOK
                % Accept as-is
                [res_supp, maxOff, idxOff] = support_debug_metrics(S, v_k, supp0);
                if ~isempty(idxOff)
                    fprintf('   → accept |supp|=%d | Sv(full)=%.2e | Sv(supp)=%.2e | max off |v|=%.2e @ %d (%s)\n', ...
                        numel(supp0), res_full, res_supp, maxOff, idxOff, rxnNames{idxOff});
                else
                    fprintf('   → accept |supp|=%d | Sv(full)=%.2e | Sv(supp)=%.2e | no reaction outside the support has flux\n', ...
                        numel(supp0), res_full, res_supp, maxOff);
                end
                
                [hasAnyBP, pairIdx] = badpairs_in_support(supp0, badPartner, bpairs);
                if hasAnyBP
                    
                    badpair_counts(pairIdx) = badpair_counts(pairIdx) + 1;  % one per pair present
                end

                if ~no_bad_pairs(supp0, badPartner)
                    efm_with_badpair_cnt = efm_with_badpair_cnt + 1;
                    continue 
                end 

                v_full = zeros(n,1);
                v_full(supp0) = v_k(supp0);
                allEFMs{i}     = supp0(:);
                allX{i}        = v_full;
                idxTwin = badPartner(supp0);
                idxTwin = idxTwin(idxTwin > 0);
                idxCover = supp0;%unique([supp0(:); idxTwin(:)]);   % twin-closure
                covered(idxCover) = true;
                found = true; break;
            end
            
            % If nullity failed, try single-drop repair on the support
            if ~nullOK || ~full_suppOK
                [ok_fix, supp_fix, v_fix] = prune_nullity_by_single_drop( ...
                    S, ub, supp0, eps_flux, tol_balance, badPartner, i);
                res_supp= norm(S(:,supp_fix) * v_fix(supp_fix), Inf);
                full_suppOK=(res_supp <= 1e-5);
                if ~full_suppOK
                    disp("sickim")
                end 
                if ok_fix && full_suppOK
                    % Accept repaired EFM
                    [res_supp, maxOff, idxOff] = support_debug_metrics(S, v_fix, supp_fix);
                    fprintf('   → accept (single-drop) |supp|=%d | Sv(full)=%.2e | Sv(supp)=%.2e | max off |v|=%.2e', ...
                            numel(supp_fix), norm(S*v_fix,Inf), res_supp, maxOff);
                    if ~isempty(idxOff), fprintf(' @ %d (%s)\n', idxOff, rxnNames{idxOff}); else, fprintf('\n'); end
                       
                    if ~no_bad_pairs(supp_fix, badPartner)
                        efm_with_badpair_cnt = efm_with_badpair_cnt + 1;
                        disp("bp")
                        continue 
                    end 
                    
                    allEFMs{i}      = supp_fix(:);
                    allX{i}         = v_fix(:);
                    idxTwin = badPartner(supp_fix);
                    idxTwin = idxTwin(idxTwin > 0);
                    idxCover = supp_fix;%unique([supp_fix(:); idxTwin(:)]);   % twin-closure
                    covered(idxCover) = true;
                    found = true; break;
                end
            end
            
            % Otherwise ban this support and try next candidate
            localTemplate = applyNoGoodCut(localTemplate, supp0, n);
            continue;

        end

        if found, break; end
    end

    if ~found
        warning('   No valid support for reaction %d (%s) after %d attempts', ...
                i, rxnNames{i}, maxSeedAttempts);
        skipped(i) = true;  % <-- avoid infinite cycling
    end
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
        break;  % exits the while-loop
    end
end

%%
missingReactions = find(~covered);
fprintf('Uncovered reactions (%d):\n', numel(missingReactions));
disp(rxnNames(missingReactions));

%% 7) Report & save — columns = EFMs, rows = [residual; reaction fluxes]
validMask = ~cellfun('isempty', allX);
efmAnchors = find(validMask);
numEFMs    = numel(efmAnchors);

if numEFMs == 0
    warning('No EFMs accepted; nothing to save.');
else
    % Build the fixed-size matrix: first row = ||S*v||_inf, then n reaction rows
    EFM_matrix = zeros(n+1, numEFMs);
    for c = 1:numEFMs
        v = allX{efmAnchors(c)}(:);        % full-length vector (n x 1), zeros where inactive
        EFM_matrix(2:end, c) = v;
        EFM_matrix(1, c)     = norm(S*v, Inf);  % steady-state residual for that EFM
    end

    % Row names = residual label + reaction names
    rowNames = [{'Sv_inf_residual'}; rxnNames(:)];

    % Column names = EFM index + (optional) anchor reaction name
    rawColNames = arrayfun(@(c) sprintf('EFM_%d_anchor_%s', c, rxnNames{efmAnchors(c)}), ...
                           1:numEFMs, 'UniformOutput', false);
    % Make them valid & unique for table VariableNames
    varNames = matlab.lang.makeValidName(rawColNames, 'ReplacementStyle','delete');
    varNames = matlab.lang.makeUniqueStrings(varNames);

    % Wrap in a table so rows "carry" the reaction names
    EFM_table = array2table(EFM_matrix, 'RowNames', rowNames, 'VariableNames', varNames);

    % Optional: binary support mask
    EFM_support = EFM_matrix(2:end, :) ~= 0;

    % Save MAT and a CSV (CSV includes row names)
    save('efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5.mat', 'EFM_matrix', 'EFM_table', 'rowNames', ...
         'varNames', 'EFM_support', 'efmAnchors', 'rxnNames', 'S', 'model_ir');

    writetable(EFM_table, 'efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5.csv', 'WriteRowNames', true);

    fprintf('Saved EFM matrix: (%d x %d) = [residual; %d reactions] x %d EFMs\n', ...
            size(EFM_matrix,1), size(EFM_matrix,2), n, numEFMs);
end

%%
load('efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5.mat','EFM_matrix','EFM_table','rowNames','varNames','EFM_support','efmAnchors','rxnNames','S','model_ir');

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
% Try removing ONE reaction (prefer not dropping 'anchor'):
%   For each j in supp0:
%     - trial = supp0 \ {j}
%     - check bad pairs
%     - LP: S(:,trial)*v=0, eps<=v<=U
%     - if feasible, re-check nullity(trial)==1 and FULL Sv<=tol
%     - accept if ok (anchor must remain)
    ok = false; supp_out = []; v_full = [];

    if numel(supp0) <= 1, return; end

    % Order: try non-anchor first, anchor last
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

    % support = thresholded OR winner-forced
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