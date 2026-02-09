% EFM_enumerator.m
% REFACTORED: Uses config.json for paths and parameters.
% with immediate no‐good cuts for invalid supports only,
% vector bad‐pair test, seeding with retries, pooling without cardinality cuts

%% 0) SETUP & CONFIG
currentScriptPath = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(currentScriptPath));
% src / efm->src->root configFile =
    fullfile(projectRoot, 'config', 'config.json');
if
  ~exist(configFile, 'file'), error('Config file not found: %s', configFile);
end config = jsondecode(fileread(configFile));
resolve = @(p) fullfile(projectRoot, strrep(p, '/', filesep));

%% 0b) Gurobi threads & paths
nThreads = feature('numcores')-2;
fprintf('Using %d threads for Gurobi\n', nThreads);
setenv('OMP_NUM_THREADS', num2str(nThreads));
if
  ~isempty(getenv('GUROBI_HOME'))
      addpath(fullfile(getenv('GUROBI_HOME'), 'examples', 'matlab'), '-begin');
end

%% 1) Load irreversible pruned model
modelPath = resolve(fullfile(config.paths.models_dir, config.model.model_file));
fprintf('Loading model: %s\n', modelPath);
data = load(modelPath);
if isfield (data, 'pruned_ir')
  model_ir = data.pruned_ir;
else
  fn = fieldnames(data);
model_ir = data.(fn{1});
end S = model_ir.S;
rxnNames = model_ir.rxns;
[ m, n ] = size(S);

% Biomass metabolite display(optional check) %
    Assuming indices from original script(2219) might be model specific,
    keeping safety check if size (S, 1) >= 2219 && S(2219, 2219) ~ =
        0 % Heuristic check met_indices = find(model_ir.S(
                                                   :, 2219) ~ = 0);
coefficients = full(model_ir.S(met_indices, 2219));
met_ids = model_ir.mets(met_indices);
if isfield (model_ir, 'metNames')
  met_names = model_ir.metNames(met_indices);
else
  met_names = met_ids;
end biomass_structure = table(met_ids, met_names, coefficients, 'VariableNames',
                              {'ID', 'Name', 'Stoichiometry'});
disp(biomass_structure);
end

%% 2) Parameters
eps_flux      = 1e-7;
tol_balance = 1e-8;
M = 5e3;
badPair_count = 0;
JACCARD_TAU = 0.99;
MAX_PER_ANCHOR = 40;
MIN_PER_ANCHOR = 3;

accepted_supps = {};
EFM_supps = {};
EFM_vecs = {};
EFM_anchor = [];

%% Compute range of all reactions for fixing the flux 
fvaModel.A        = sparse(S);
fvaModel.rhs = zeros(m, 1);
fvaModel.sense = repmat('=', m, 1);
fvaModel.lb = model_ir.lb;
fvaModel.ub = model_ir.ub;
fvaModel.vtype = repmat('C', n, 1);
fvaModel.modelsense = 'min';
fvaModel.obj = zeros(n, 1);

% Gurobi parameters for FVA
fvaParams.OutputFlag = 0;
fvaParams.TimeLimit = 60;
fvaParams.FeasibilityTol = 1e-9;

vMin = nan(n, 1);
vMax = nan(n, 1);

fprintf('Running FVA Pre-check...\n');
for j = 1:n
    % 1) minimize v_j
    fvaModel.obj(:) = 0;
fvaModel.obj(j) = 1;
fvaModel.modelsense = 'min';
sol = gurobi(fvaModel, fvaParams);
if isfield (sol, 'status')
  &&strcmp(sol.status, 'OPTIMAL'), vMin(j) = sol.objval;
else
  , warning('FVA min infeasible %d', j); end

    % 2) maximize v_j
    fvaModel.modelsense = 'max';
sol = gurobi(fvaModel, fvaParams);
if isfield (sol, 'status')
  &&strcmp(sol.status, 'OPTIMAL'), vMax(j) = sol.objval;
else
  , warning('FVA max infeasible %d', j);
end end range = [ vMin, vMax ];

idxZeroMax = abs(range( :, 2)) <= 1e-5;
fprintf('Reactions with vMax <= 1e-5: %d\n', nnz(idxZeroMax));

%% 3) Build static MILP template 
vMin(~isfinite(vMin)) = 0;
vMax(~isfinite(vMax)) = M;
vMin = max(vMin, 0);
range = [ vMin, vMax ];

lb = [range( :, 1); zeros(n, 1)];
ub = [range( :, 2); ones(n, 1)];
ub_eff = ub;
ub_eff(~isfinite(ub_eff)) = M;
ub = ub_eff;

vtype = [repmat('C', n, 1); repmat('B', n, 1)];

Aeq = [ sparse(S), sparse(m, n) ];
beq = zeros(m, 1);

A1 = [ speye(n), -ub(1 : n).*speye(n) - speye(n), eps_flux *speye(n) ];
b1 = zeros(2 * n, 1);

nz = (S ~ = 0);
importCols = false(1, n);
exportCols = false(1, n);
for
  j = 1 : n rows = nz( :, j);
if any (rows)
  col = S(rows, j);
importCols(j) = all(col > 0);
exportCols(j) = all(col < 0);
end end Aimp = [];
bimp = [];
if any (importCols)
  , r = sparse(1, 2 * n);
r(n + find(importCols)) = -1;
Aimp = r;
bimp = -1;
end

    Aexp = [];
bexp = [];
if any (exportCols)
  , r = sparse(1, 2 * n);
r(n + find(exportCols)) = -1;
Aexp = r;
bexp = -1;
end

    % Bad Pairs(Reversible split handling) badPartner = zeros(n, 1);
badPair_count = 0;
for
  i = 1 : n rn = rxnNames{i};
if endsWith (rn, '_f')
  j = find(strcmp(rxnNames, [rn(1 : end - 2) '_b']), 1);
if
  ~isempty(j) badPartner(i) = j;
badPartner(j) = i;
if i
  < j, badPair_count = badPair_count + 1; end
        end
    end
end
bpairs = find(badPartner>0 & (1:n)' < badPartner);
p = numel(bpairs);
Abp = sparse(p,n*2);
bbp = ones(p,1);
for k=1:p
    i = bpairs(k); j = badPartner(i);
    Abp(k, n+i) = 1; Abp(k, n+j) = 1;
end

Aineq = [ A1; Aimp; Aexp; Abp];
bineq = [ b1; bimp; bexp; bbp];
sense = [ repmat('=', size(Aeq,1),1); repmat('<', size(Aineq,1),1) ];
obj = [ zeros(n,1) ; ones(n,1) ];

linModel.A = [Aeq; Aineq]; linModel.rhs = [beq; bineq]; linModel.sense = sense;
linModel.lb = lb; linModel.ub = ub; linModel.vtype = vtype; linModel.modelsense = 'min'; linModel.obj = obj;

%% 4) Gurobi parameters
commonParams = struct('FeasibilityTol',1e-9, 'IntFeasTol', 1e-9, 'NumericFocus', 3);
pPool = struct('OutputFlag', 0, 'PoolSearchMode', 2, 'PoolSolutions', 100,
               'TimeLimit', 60, 'MIPGap', 0.0, 'Cuts', 1, 'Presolve', 2,
               'Heuristics', 0.5, 'MIPFocus', 1, 'Threads', nThreads);
pSeed = struct('OutputFlag', 0, 'PoolSearchMode', 2, 'PoolSolutions', 40,
               'TimeLimit', 600, 'MIPGap', 0.3, 'Cuts', -1, 'Presolve', 2,
               'Heuristics', 0.2, 'MIPFocus', 1, 'Threads', nThreads);
for
  p_ = {pPool, pSeed}, p_{1}.FeasibilityTol = commonParams.FeasibilityTol;
p_{1}.IntFeasTol = commonParams.IntFeasTol;
p_{1}.NumericFocus = commonParams.NumericFocus;
end

    fvaParams.OutputFlag = 0;
fvaParams.TimeLimit = 60;
fvaParams.FeasibilityTol = 1e-9;

%% 5) Checkpoint / Resume (timed saves only)
% Save dir logic
efmOutDir = fullfile(resolve(config.paths.models_dir), 'E_coli', 'efms');
% Assuming structure if ~exist(efmOutDir, 'dir'), mkdir(efmOutDir);
end

    CHECKPOINT_ENABLE = true;
AUTOSAVE_SEC = 15 * 60;
CHECKPOINT_FILE = fullfile(efmOutDir, 'efm_checkpoint.mat');

LOAD_CHECKPOINT = false;
LOAD_FILE = CHECKPOINT_FILE;
model_sig = struct('m', m, 'n', n, 'nnzS', nnz(S), 'rxnCount', numel(rxnNames));

if LOAD_CHECKPOINT
  &&exist(LOAD_FILE, 'file') R = load(LOAD_FILE, 'model_sig');
if isfield (R, 'model_sig')
  &&isequal(R.model_sig, model_sig)
      fprintf('[resume] Loading: %s\n', LOAD_FILE);
load_checkpoint_state(LOAD_FILE);
if exist ('rng_state', 'var')
  , rng(rng_state);
end else warning('[resume] Incompatible checkpoint; starting fresh.');
end end

    if ~exist('covered', 'var'),
    covered = false(n, 1);
end if ~exist('skipped', 'var'), skipped = false(n, 1);
end

    lastSave = tic;
targetBiomassName = 'BIOMASS_Ec_iML1515_WT_75p37M';
bio_idx = find(strcmp(model_ir.rxns, targetBiomassName));

%% 6) Enumeration loop
covered  = false(n,1);
skipped = false(n, 1);
allEFMs = {};
allX = {};
prevUc = 0;
efm_count_total = 0;
efm_with_badpair_cnt = 0;
badpair_counts = zeros(p, 1);

maxSeedAttempts = 25;
stallLimit = 300;
noProgressCount = 0;
prevCoveredCount = sum(covered);
stallTerminated = false;

pool = gcp('nocreate');
if isempty (pool)
  , pool = parpool('local', nThreads);
end

        try while any(~covered & ~skipped) if CHECKPOINT_ENABLE &&toc(
            lastSave) > AUTOSAVE_SEC
    save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
lastSave = tic;
end

    toTry = find(~covered & ~skipped);
fprintf('\n=== Next reaction: %d remaining (%d covered) ===\n', numel(toTry),
        sum(covered));
idx = randi(numel(toTry));
    % Priority logic (Biomass first if uncovered?)
    if ~covered(bio_idx) && ismember(bio_idx, toTry), idx = find(toTry==bio_idx);
    end

        i = toTry(idx);
    fprintf('I: %d (%s)\n', i, rxnNames{i});
    prevUc = i;

    localTemplate = linModel;
    found = false;
    accepted_this_anchor = 0;

    for
      attempt = 1 : maxSeedAttempts M1 = localTemplate;
    M1.lb(n + i) = 1;
    M1.ub(n + i) = 1;
    M1.lb(i) = 1;
    sol1 = gurobi(M1, pSeed);
    if
      ~isfield(sol1, 'pool') || isempty(sol1.pool), continue;
    end

        numCand = numel(sol1.pool);
    fprintf('  attempt %d: %d pool candidates\n', attempt, numCand);

    % Pre - allocate acceptType = zeros(numCand, 1);
    supp_accept = cell(numCand, 1);
    v_accept = cell(numCand, 1);
    nogood_supp = cell(numCand, 1);
    anchorInSupp = false(numCand, 1);
    hasPair_list = false(numCand, 1);
    pairIdx_list = cell(numCand, 1);
    badpair_violate = false(numCand, 1);

    parfor k = 1 : numCand xk = getSolutionVector(sol1.pool(k));
    v_k = xk(1 : n);
    supp0 = find(v_k >= eps_flux - tol_balance);
    if isempty (supp0)
      || ~ismember(i, supp0), continue;
    end if ~no_bad_pairs(supp0, badPartner)[supp0, v_k] =
        cancel_twins(supp0, v_k, badPartner, S, eps_flux, tol_balance);
    if isempty (supp0)
      || ~ismember(i, supp0), continue;
    end end anchorInSupp(k) = true;
    nogood_supp{k} = supp0;
    [ hBP, pIdx ] = badpairs_in_support(supp0, badPartner, bpairs);
    hasPair_list(k) = hBP;
    if hBP
      , pairIdx_list{k} = pIdx;
    end

        nullOK = nullity1_ok(S, supp0);
    res_supp = norm(S( :, supp0) * v_k(supp0), Inf);
    full_suppOK = (res_supp <= 1e-8);
    stillBP = ~no_bad_pairs(supp0, badPartner);
    badpair_violate(k) = stillBP;

    if nullOK
      &&full_suppOK && ~stillBP acceptType(k) = 1;
    supp_accept{k} = supp0;
    v_accept{k} = v_k;
    continue;
    end

        if (~nullOK || ~full_suppOK)[ok_fix, supp_fix, v_fix] =
            prune_nullity_by_single_drop(S, ub, supp0, eps_flux, tol_balance,
                                         badPartner, i);
    if ok_fix
      res_supp_fix = norm(S( :, supp_fix) * v_fix(supp_fix), Inf);
    full_suppOK_fix = (res_supp_fix <= 1e-8);
    stillBP_fix = ~no_bad_pairs(supp_fix, badPartner);
    if full_suppOK_fix
      &&~stillBP_fix acceptType(k) = 2;
    supp_accept{k} = supp_fix( :);
    v_accept{k} = v_fix( :);
                    end
                end
            end
        end

        for k = 1:numCand
            if ~anchorInSupp(k), continue;
                    end if hasPair_list (k) && ~isempty(pairIdx_list{k}),
                        badpair_counts(pairIdx_list{
                            k}) = badpair_counts(pairIdx_list{k}) + 1;
                    end if badpair_violate (k),
                        efm_with_badpair_cnt = efm_with_badpair_cnt + 1;
                    continue;
                    end

                            if acceptType (k) == 1 ||
                        acceptType(k) == 2 supp0 = supp_accept{k};
                    v_k = v_accept{k};
                    tag = '';
                    if acceptType (k)
                      == 2, tag = ' (single-drop)';
                    end

                        J = 0;
                    if
                      ~isempty(accepted_supps), J = jaccard_max_sim(
                                                    supp0, accepted_supps);
                    end if J >= JACCARD_TAU,
                        localTemplate = applyNoGoodCut(localTemplate, supp0, n);
                    continue;
                    end

                        [okLP, v_full2, ~] =
                            maximize_sum_v_on_support(S, ub, supp0, eps_flux);
                    if okLP
                      , v_full = v_full2;
                    else
                      , v_full = zeros(n, 1);
                    v_full(supp0) = v_k(supp0);
                    end

                        EFM_supps{end + 1} = supp0(
                            :);
                    EFM_vecs{end + 1} = v_full( :);
                    EFM_anchor(end + 1) = i;
                    accepted_supps{end + 1} = supp0( :);
                    covered(supp0) = true;
                    localTemplate = applyNoGoodCut(localTemplate, supp0, n);
                    fprintf(
                        '   → accept%s |supp|=%d | Sv=%.2e | sum v=%.3g | J=%.2f\n',
                        tag, numel(supp0), norm(S *v_full, Inf),
                        sum(v_full(supp0)), J);
                    accepted_this_anchor = accepted_this_anchor + 1;
                    localTemplate =
                        applyNoGoodCut(localTemplate, nogood_supp{k}, n);

                    if CHECKPOINT_ENABLE
                      &&toc(lastSave) > AUTOSAVE_SEC,
                          save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
                    lastSave = tic;
                    end else localTemplate =
                        applyNoGoodCut(localTemplate, nogood_supp{k}, n);
                    end end

                            if accepted_this_anchor >= MAX_PER_ANCHOR ||
                        accepted_this_anchor >= MIN_PER_ANCHOR,
                        break;
                    end end

                        if accepted_this_anchor ==
                        0 warning('   No valid support logic...');
                    skipped(i) = true;
                    end

                        currCoveredCount = sum(covered);
                    if currCoveredCount
                      > prevCoveredCount, noProgressCount = 0;
                    else
                      , noProgressCount = noProgressCount + 1;
                    end prevCoveredCount = currCoveredCount;

                    if noProgressCount
                      >= stallLimit fprintf(
                             '\nNo new coverage for %d consecutive iterations. Early termination.\n',
                             stallLimit);
                    stallTerminated = true;
                    if CHECKPOINT_ENABLE
                      , save_checkpoint_atomic(CHECKPOINT_FILE, model_sig);
                    lastSave = tic;
                    end break;
                    end end catch ME fprintf('[error] %s\n', ME.message);
                    rethrow(ME);
end

%% 7) Report & save
numEFMs = numel(EFM_vecs);
if numEFMs
  > 0 EFM_matrix = zeros(n + 1, numEFMs);
    for
      c = 1 : numEFMs v = EFM_vecs{c}( :);
    EFM_matrix(2 : end, c) = v;
    EFM_matrix(1, c) = norm(S * v, Inf);
    end rowNames = [ {'Sv_inf_residual'}; rxnNames( :) ];
    rawColNames =
        arrayfun(@(c) sprintf('EFM_%d_anchor_%s', c, rxnNames{EFM_anchor(c)}), 1
                 : numEFMs, 'UniformOutput', false);
    varNames =
        matlab.lang.makeValidName(rawColNames, 'ReplacementStyle', 'delete');
    varNames = matlab.lang.makeUniqueStrings(varNames);
    EFM_table = array2table(EFM_matrix, 'RowNames', rowNames, 'VariableNames',
                            varNames);
    EFM_support = abs(EFM_matrix(2 : end, :)) > (eps_flux - tol_balance);

    saveFile = fullfile(efmOutDir, 'efms_matrix_iML1515_forBiomass2.mat');
    save(saveFile, 'EFM_matrix', 'EFM_table', 'rowNames', 'varNames',
         'EFM_support', 'EFM_anchor', 'EFM_supps', 'rxnNames', 'S', 'model_ir');

    csvFile = fullfile(efmOutDir, 'efms_matrix_iML1515_forBiomass2.csv');
    writetable(EFM_table, csvFile, 'WriteRowNames', true);
    fprintf('Saved: %s\n', saveFile);
    end

        % % — Helpers(unchanged except context) — function x =
        getSolutionVector(sol) if isfield (sol, 'x') && ~isempty(sol.x),
                                                           x = sol.x;
    elseif isfield(sol, 'xn') && ~isempty(sol.xn), x = sol.xn;
    elseif isfield(sol, 'pool') && ~isempty(sol.pool) &&
        isfield(sol.pool(1), 'xn'),
        x = sol.pool(1).xn;
    else, error('No solution vector found.');
    end end function T = applyNoGoodCut(T, supp, n) oth = setdiff(1 : n, supp);
    cut = sparse(1, 2 * n);
    cut(n + supp) = 1;
    cut(n + oth) = -1;
    T.A = [T.A; cut];
    T.rhs = [T.rhs; numel(supp) - 1];
    T.sense = [T.sense; '<'];
    end function tf = no_bad_pairs(supp, badPartner), tf = true;
    if isempty (supp)
      , return;
    end;
    bp = badPartner(supp);
    if any (bp > 0)
      , tf = ~any(ismember(bp(bp > 0), supp));
    end;
    end function supp_out =
        simplify_bad_pairs(supp_in, badPartner, keepIdx) supp = supp_in(
            :)'; changed = true; while changed changed = false;
        for
          t = 1 : numel(supp) i = supp(t);
        j = badPartner(i);
        if j
          > 0 && ismember(j, supp) if i == keepIdx, drop = j;
        elseif j == keepIdx, drop = i;
        else, drop = max(i, j);
        end supp(supp == drop) = [];
        changed = true;
        break;
        end end end supp_out = supp( :);
        end function ok = nullity1_ok(S, supp), if isempty (supp), ok = false;
        return;
        end;
        [ ~, R ] = qr(S( :, supp), 0);
        rk = sum(abs(diag(R)) > 1e-10);
        ok = (numel(supp) - rk == 1);
        end function[ok, supp_out, v_full] = prune_nullity_by_single_drop(
            S, ub_eff, supp0, eps_flux, tol_balance, badPartner, anchor) ok =
            false;
        supp_out = [];
        v_full = [];
        if numel (supp0)
          <= 1, return; end
    order = supp0(:).'; if nargin >= 7 && ~isempty(anchor), order = [ order(order ~= anchor), anchor ]; end
    for j = order
        trial = setdiff(supp0, j, 'stable');
        if isempty (trial)
          , continue;
        end if ~no_bad_pairs(trial, badPartner), continue;
        end if nargin >= 7 && ~isempty(anchor) && ~ismember(anchor, trial),
            continue;
        end[ok_refit, v_trial] =
            refit_on_support_feasible(S, ub_eff, trial, eps_flux);
        if
          ~ok_refit, continue;
        end if ~nullity1_ok(S, trial), continue;
        end n = size(S, 2);
        v_full = zeros(n, 1);
        v_full(trial) = v_trial;
        if norm (S( :, trial) * v_full(trial), Inf)
          > 1e-5, continue;
        end ok = true;
        supp_out = trial( :);
        return;
        end end function[ok, v_s] =
            refit_on_support_feasible(S, ub_eff, supp, eps_pos) m = size(S, 1);
        k = numel(supp);
        if k
          == 0, ok = false;
        v_s = [];
        return;
        end model.A = sparse(S( :, supp));
        model.rhs = zeros(m, 1);
        model.sense = repmat('=', m, 1);
        model.lb = eps_pos * ones(k, 1);
        model.ub = ub_eff(supp);
        model.vtype = repmat('C', k, 1);
        model.modelsense = 'min';
        model.obj = zeros(k, 1);
        sol = gurobi(model, struct('OutputFlag', 0, 'FeasibilityTol', 1e-9));
        ok = isfield(sol, 'status') && strcmp(sol.status, 'OPTIMAL');
        if ok
          , v_s = sol.x;
        else
          , v_s = [];
        end end function[res_supp, maxOff, idxOff] =
            support_debug_metrics(S, v_raw, supp) if isempty (supp),
                                           res_supp = 0;
        else, res_supp = norm(S( :, supp) * v_raw(supp), Inf);
        end n = size(S, 2);
        off = setdiff(1 : n, supp);
        if isempty (off)
          , maxOff = 0;
        idxOff = [];
        return;
        end[maxOff, k] = max(abs(v_raw(off)));
        idxOff = off(k);
        end function[supp_out, v_new] =
            cancel_twins(supp_in, v, badPartner, S, eps, tol) v_new = v(
                :); supp  = supp_in(:).'; n = numel(v_new); keepWinner = false(n,1); seen = false(n,1);
    for t = 1:numel(supp)
        i = supp(t);
        j = badPartner(i);
        if j
          <= 0 || ~ismember(j, supp) || seen(i) || seen(j), continue;
        end seen([i j]) = true;
        if max (abs(S( :, i) + S( :, j)))
          <= 1e-9 vi = v_new(i);
        vj = v_new(j);
        if vi
          > 0 && vj > 0 if vi == vj, keepWinner(i) = true;
        keepWinner(j) = true;
        continue;
        end d = min(vi, vj);
        v_new(i) = vi - d;
        v_new(j) = vj - d;
        if v_new (i)
          > v_new(j) && v_new(i) > tol, keepWinner(i) = true;
        elseif v_new(j) > tol, keepWinner(j) = true;
        end end end end supp_out =
            unique([find(v_new >= eps - tol); find(keepWinner)]);
        end function[hasAny, idxPairs] =
            badpairs_in_support(supp, badPartner, bpairs) if isempty (supp),
                             hasAny = false;
        idxPairs = [];
        return;
        end inSupp = false(size(badPartner));
        inSupp(supp) = true;
        j = badPartner(bpairs);
        present = inSupp(bpairs) & inSupp(j);
        hasAny = any(present);
        idxPairs = find(present);
        end function[ok, v_full, objval] =
            maximize_sum_v_on_support(S, ub_eff, supp, eps_pos) m = size(S, 1);
        n = size(S, 2);
        k = numel(supp);
        if k
          == 0, ok = false;
        v_full = [];
        objval = -Inf;
        return;
        end model.A = sparse(S( :, supp));
        model.rhs = zeros(m, 1);
        model.sense = repmat('=', m, 1);
        model.lb = eps_pos * ones(k, 1);
        model.ub = ub_eff(supp);
        model.vtype = repmat('C', k, 1);
        model.modelsense = 'max';
        model.obj = ones(k, 1);
        sol = gurobi(model, struct('OutputFlag', 0, 'FeasibilityTol', 1e-9));
        ok = isfield(sol, 'status') && strcmp(sol.status, 'OPTIMAL');
        if ok
          , v_full = zeros(n, 1);
        v_full(supp) = sol.x;
        objval = sol.objval;
        else, v_full = [];
        objval = -Inf;
        end end function J = jaccard_max_sim(supp, accepted_supps) J = 0;
        A = unique(supp( :));
    for
      t = 1 : numel(accepted_supps) B = unique(accepted_supps{t}( :));
    inter = numel(intersect(A, B));
    uni = numel(union(A, B));
    j = inter / max(uni, 1);
    if j
      > J, J = j;
    end if J >= 0.999, return;
    end end end function save_checkpoint_atomic(targetFile, model_sig)
    try rng_state = rng;
    catch, rng_state = [];
    end varNames = {'EFM_supps',
                    'EFM_vecs',
                    'EFM_anchor',
                    'accepted_supps',
                    'covered',
                    'skipped',
                    'linModel',
                    'noProgressCount',
                    'prevCoveredCount',
                    'stallTerminated',
                    'efm_with_badpair_cnt',
                    'badpair_counts',
                    'JACCARD_TAU',
                    'MAX_PER_ANCHOR',
                    'MIN_PER_ANCHOR',
                    'eps_flux',
                    'tol_balance',
                    'M',
                    'ub',
                    'badPartner',
                    'bpairs',
                    'p',
                    'S',
                    'rxnNames',
                    'model_sig',
                    'rng_state'};
    Sstruct = struct();
    for
      v = varNames, vn = v{1};
    if evalin ('caller', sprintf('exist(' '%s' ',' 'var' ')', vn))
      , Sstruct.(vn) = evalin('caller', vn);
    end;
    end tmp = [ targetFile, '.tmp' ];
    save(tmp, '-struct', 'Sstruct', '-v7.3');
    movefile(tmp, targetFile, 'f');
    fprintf('[checkpoint] %s\n', targetFile);
    end function load_checkpoint_state(fname) L = load(fname);
    f = fieldnames(L); for
      k = 1 : numel(f), assignin('caller', f{k}, L.(f{k}));
    end fprintf('[resume] state loaded from %s\n', fname);
    end
