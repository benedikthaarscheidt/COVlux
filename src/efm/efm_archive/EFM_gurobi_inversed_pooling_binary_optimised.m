% with immediate no‐good cuts for invalid supports only,
% vector bad‐pair test, seeding with retries, pooling without cardinality cuts

%% 0) Gurobi threads & paths
nThreads = feature('numcores');
fprintf('Using %d threads for Gurobi\n',nThreads);
setenv('OMP_NUM_THREADS',num2str(nThreads));
addpath(fullfile(getenv('GUROBI_HOME'),'examples','matlab'),'-begin');

%% 1) Load irreversible pruned model
data     = load('/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/e_coli_core_splitandpruned.mat','pruned_ir')
model_ir = data.pruned_ir;
S        = model_ir.S;
rxnNames = model_ir.rxns;
[m,n]    = size(S);

%% 2) Parameters
eps_flux      = 1e-5;
tol_balance   = 1e-6;
M             = 1e3;
badPair_count = 0;

%% 3) Build static MILP template (bare template)
aeq = [S, sparse(m,n)];
beq = zeros(m,1);
A1  = [ eye(n), -M*eye(n)
       -eye(n),  eps_flux*eye(n) ];
b1  = zeros(2*n,1);

% import/export cuts
nz         = (S~=0);
importCols = false(1,n);
exportCols = false(1,n);
for j = 1:n
    rows = nz(:,j);
    if any(rows)
        importCols(j) = all(S(rows,j)>0);
        exportCols(j) = all(S(rows,j)<0);
    end
end
cutImp = sparse(1,2*n); cutImp(n+find(importCols)) = -1;
cutExp = sparse(1,2*n); cutExp(n+find(exportCols)) = -1;

Aineq = [A1; cutImp; cutExp];
bineq = [b1; -1; -1];

lb0    = [model_ir.lb; zeros(n,1)];
ub0    = [model_ir.ub; ones(n,1)];
vtype0 = [repmat('C',n,1); repmat('B',n,1)];

bareTemplate.A          = [aeq; Aineq];
bareTemplate.rhs        = [beq; bineq];
bareTemplate.sense      = [ repmat('=',m,1); repmat('<',size(Aineq,1),1) ];
bareTemplate.lb         = lb0;
bareTemplate.ub         = ub0;
bareTemplate.vtype      = vtype0;
bareTemplate.modelsense = 'min';
bareTemplate.obj        = [ zeros(n,1) ; ones(n,1) ];

%% 4) Gurobi parameters
commonParams = struct( ...
   'FeasibilityTol',1e-9, ...
   'IntFeasTol',    1e-7, ...
   'NumericFocus',  1);
pPool = struct( ...
  'OutputFlag',     0, ...
  'PoolSearchMode',2, ...
  'PoolSolutions',1, ...
  'TimeLimit',    60, ...
  'MIPGap',0.5, ...
  'Cuts',1, ...
  'Presolve',2, ...%this was 2 
  'Heuristics',0.5, ...
  'MIPFocus',1, ...
  'Threads',nThreads );
pSeed = struct( ...
  'OutputFlag',     0, ...
  'PoolSearchMode',0, ...
  'PoolSolutions',1, ...
  'TimeLimit',    600, ...
  'MIPGap',0.01, ...
  'Cuts',1, ... 
  'Presolve',2, ...
  'Heuristics',0.2, ...
  'MIPFocus',1, ...
  'Threads',nThreads );
for p = {pPool,pSeed}
    p{1}.FeasibilityTol = commonParams.FeasibilityTol;
    p{1}.IntFeasTol     = commonParams.IntFeasTol;
    p{1}.NumericFocus   = commonParams.NumericFocus;
end

%% 5) Precompute bad‐pair map
badPartner = zeros(n,1);
for i = 1:n
    rn = rxnNames{i};
    if endsWith(rn,'_f')
        j = find(strcmp(rxnNames,[rn(1:end-2) '_b']),1);
        if ~isempty(j)
            badPartner(i)=j;
            badPartner(j)=i;
        end
    end
end

%% 6) Enumeration loop
covered  = false(n,1);
skipped  = false(n,1);
allEFMs  = {};
allX     = {};
prevUc   = 0;

failCounts = struct('nullity',0,'balance',0,'fluxpos',0,'badpair',0,'anchor',0,'emptysupp',0);
failSizes  = struct('nullity',[],'balance',[],'fluxpos',[],'badpair',[],'anchor',[],'emptysupp',[]);




maxSeedAttempts = 1;
% start parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local',nThreads);
end

% initialize globalTemplate (starts bare, then accumulates only valid EFM cuts)
globalTemplate = bareTemplate;

while any(~covered & ~skipped ) %
    batch = sum(skipped | ~covered);
    toTry = find(~covered & ~skipped);
    fprintf('\n=== Next reaction: %d remaining (%d covered) ===\n', numel(toTry), sum(covered));

    % pick next uncovered reaction
    idx = find(toTry > prevUc,1,'first');
    if isempty(idx), idx = 1; end
    uc = toTry(idx);
    prevUc = uc;

    % SEEDING: use a fresh local template per reaction
    localTemplate = globalTemplate;
    acceptedSeeds = 0;

    for attempt = 1:maxSeedAttempts
        % build & solve seeded MILP
        M1 = localTemplate;
        cutF = sparse(1,2*n); cutF(n+uc) = 1;
        M1.A      = [M1.A;   cutF];
        M1.rhs    = [M1.rhs; 1];
        M1.sense  = [M1.sense; '>'];
        M1.lb(n+uc) = 1;
        sol1 = gurobi(M1, pSeed);

        % extract pool seeds
        
    if isfield(sol1,'pool') && ~isempty(sol1.pool)
            poolSeeds = sol1.pool;
        elseif isfield(sol1,'xn')
            poolSeeds = struct('xn',sol1.xn);
        else
            poolSeeds = struct('xn',sol1.x);
        end

        % parallel nullity‐check & collect supports
        numCand  = numel(poolSeeds);
        valid    = false(numCand,1);
        suppList = cell(numCand,1);
        vList    = cell(numCand,1);
        parfor k = 1:numCand
            xk = getSolutionVector(poolSeeds(k));
            [ok, supp, fl] = checkSupport(xk, S, tol_balance, badPartner, eps_flux, uc);
            suppList{k} = supp;
            flList{k}   = fl;
            if ok && ismember(uc, supp)
                valid(k) = true;
                %vList{k} = xk(1:n);
                vList{k} = xk(supp);
            end
        end

        % ban only invalid supports on localTemplate
        for k = 1:numCand
            if ~valid(k)
                fl = flList{k};
                if fl.emptysupp, failCounts.emptysupp = failCounts.emptysupp+1; failSizes.emptysupp(end+1,1)=fl.supp_size; end
                if ~fl.nullity_ok, failCounts.nullity = failCounts.nullity+1; failSizes.nullity(end+1,1)=fl.supp_size; end
                if ~fl.balance_ok, failCounts.balance = failCounts.balance+1; failSizes.balance(end+1,1)=fl.supp_size; end
                if ~fl.fluxpos_ok, failCounts.fluxpos = failCounts.fluxpos+1; failSizes.fluxpos(end+1,1)=fl.supp_size; end
                if ~fl.badpair_ok, failCounts.badpair = failCounts.badpair+1; failSizes.badpair(end+1,1)=fl.supp_size; end
                if ~fl.anchor_ok,  failCounts.anchor  = failCounts.anchor +1; failSizes.anchor (end+1,1)=fl.supp_size; end
                localTemplate = applyNoGoodCut(localTemplate, suppList{k}, n);
            end
        end

        % accept valid seeds: mark covered, record and add permanent cuts
        if any(valid)
            acceptedSeeds = sum(valid);
            for k = find(valid).'
                supp = suppList{k};
                covered(supp)       = true;
                allEFMs{end+1}      = supp;
                allX{end+1}         = vList{k};
                %localTemplate      = applyNoGoodCut(globalTemplate, supp, n);
            end
            fprintf('  Seeding: accepted %d on attempt %d\n', acceptedSeeds, attempt);
            break
        else
            fprintf('  Seeding: no valid seeds on attempt %d for reaction %d\n', attempt, uc);
        end
    end

    if acceptedSeeds == 0
        warning('  Skipping reaction %d after %d seeding attempts\n', uc, maxSeedAttempts);
        skipped(uc) = true;
    end

    %% POOLING: restart from globalTemplate (only permanent cuts)
    %%localTemplate = globalTemplate;
    %M2 = localTemplate;
    %rem = find(~covered);
    %cutU = sparse(1,2*n); cutU(n+rem) = 1;
    %M2.A      = [M2.A;   cutU];
    %M2.rhs    = [M2.rhs; 1];
    %M2.sense  = [M2.sense; '>'];
    %M2.lb(n+uc) = 1;
    %sol2 = gurobi(M2, pPool);
%
    % if isfield(sol2,'pool') && ~isempty(sol2.pool)
    %    poolSol = sol2.pool;
    %elseif isfield(sol2,'xn')
    %    poolSol = struct('xn',sol2.xn);
    %else
    %    poolSol = struct('xn',sol2.x);
    %end
%
    %numPools = numel(poolSol);
    %validP    = false(numPools,1);
    %suppListP = cell(numPools,1);
    %vListP    = cell(numPools,1);
    %flListP   = cell(numPools,1);
    %parfor p = 1:numPools
    %    xk = getSolutionVector(poolSol(p));
    %    [ok, supp, fl] = checkSupport(xk, S, tol_balance, badPartner, eps_flux, []);
    %    suppListP{p} = supp;
    %    flListP{p}   = fl;
    %    if ok
    %        validP(p) = true;
    %        vListP{p} = xk(1:n);
    %    end
    %end
%
    %for p = 1:numPools
    %    if validP(p)
    %        supp_p = suppListP{p};      % numeric index vector
    %        v_p    = vListP{p};         % flux vector (n×1)
    %
    %        covered(supp_p(:)) = true;  % mark covered
    %        allEFMs{end+1}     = supp_p(:);
    %        allX{end+1}        = v_p(:);
    %    else
    %        fl = flListP{p};
    %        if fl.emptysupp, failCounts.emptysupp = failCounts.emptysupp+1; failSizes.emptysupp(end+1,1)=fl.supp_size; end
    %        if ~fl.nullity_ok, failCounts.nullity = failCounts.nullity+1; failSizes.nullity(end+1,1)=fl.supp_size; end
    %        if ~fl.balance_ok, failCounts.balance = failCounts.balance+1; failSizes.balance(end+1,1)=fl.supp_size; end
    %        if ~fl.fluxpos_ok, failCounts.fluxpos = failCounts.fluxpos+1; failSizes.fluxpos(end+1,1)=fl.supp_size; end
    %        if ~fl.badpair_ok, failCounts.badpair = failCounts.badpair+1; failSizes.badpair(end+1,1)=fl.supp_size; end
    %        if ~fl.anchor_ok,  failCounts.anchor  = failCounts.anchor +1; failSizes.anchor (end+1,1)=fl.supp_size; end
    %    end
    %end
    %fprintf('  Pooling: accepted %d modes\n', sum(valid));
end


null_sz = failSizes.nullity;
valid_sz = cellfun(@numel, allEFMs);

if isempty(null_sz)
    fprintf('\nNullity failed 0 times.\n');
else
    [sz_null,~,icn] = unique(null_sz);
    cnt_null = accumarray(icn,1);
    [sz_null,ordn] = sort(sz_null); cnt_null = cnt_null(ordn);

    [sz_val,~,icv] = unique(valid_sz);
    cnt_val = accumarray(icv,1);
    [sz_val,ordv] = sort(sz_val); cnt_val = cnt_val(ordv);

    sz_all = union(sz_null, sz_val);
    cnt_null_all = zeros(size(sz_all));
    cnt_val_all  = zeros(size(sz_all));
    [~,ia] = ismember(sz_null, sz_all); cnt_null_all(ia) = cnt_null;
    [~,ib] = ismember(sz_val,  sz_all); cnt_val_all(ib)  = cnt_val;

    rate_null_vs_valid = cnt_null_all ./ max(1, cnt_null_all + cnt_val_all);

    fprintf('\nNullity failures by support size:\n');
    for t = 1:numel(sz_null)
        fprintf('  size %d : %d\n', sz_null(t), cnt_null(t));
    end

    mn_null = mean(null_sz); md_null = median(null_sz);
    mn_val  = mean(valid_sz); md_val  = median(valid_sz);
    fprintf('\nSizes — mean/median:\n');
    fprintf('  nullity fails : %.2f / %.2f\n', mn_null, md_null);
    fprintf('  accepted EFMs : %.2f / %.2f\n', mn_val,  md_val);

    fprintf('\nNullity-fail share vs accepted by size:\n');
    for t = 1:numel(sz_all)
        fprintf('  size %d : %d null / %d valid  →  %.2f\n', ...
            sz_all(t), cnt_null_all(t), cnt_val_all(t), rate_null_vs_valid(t));
    end
end

%% 7) Report & save
uncov = find(~covered);
if ~isempty(uncov)
    fprintf('\nUncovered (%d): %s\n', numel(uncov), mat2str(uncov'));
else
    fprintf('\nAll covered.\n');
end
skList = find(skipped);
if ~isempty(skList)
    fprintf('Skipped (%d): %s\n', numel(skList), mat2str(skList'));
end
fprintf('\nFailure counts:\n');
fprintf('  nullity     : %d\n', failCounts.nullity);
fprintf('  steady-state: %d\n', failCounts.balance);
fprintf('  flux<eps    : %d\n', failCounts.fluxpos);
fprintf('  bad-pair    : %d\n', failCounts.badpair);
fprintf('  no-anchor   : %d\n', failCounts.anchor);
fprintf('  empty-supp  : %d\n', failCounts.emptysupp);

fluxMat = zeros(n,numel(allEFMs));
for i = 1:numel(allEFMs)
    v = allX{i};
    if numel(v)==n
        fluxMat(:,i) = v;
    else
        tmp = zeros(n,1); tmp(allEFMs{i}) = v; fluxMat(:,i) = tmp;
    end
end
save('EFMs_pool_badpair_nocard_parpool_flags.mat','allEFMs','fluxMat','covered','skipped','badPair_count','failCounts','failSizes');
disp("done");

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

function [isValid, supp, fl] = checkSupport(x, S, tol_balance, badPartner, eps_flux, anchorIdx)
    n      = size(S,2);
    v      = x(1:n);
    b      = x(n+1:2*n); %#ok<NASGU>
    supp   = find(v >= eps_flux);

    % --- report largest flux outside support ---
    inactive = setdiff(1:n, supp);
    if isempty(inactive)
        fl.max_inactive_flux = 0;
        fl.max_inactive_rxn  = NaN;
        fprintf('  no inactive reactions (support covers all)\n');
    else
        [max_flux, idxMax]   = max(abs(v(inactive)));
        rxn_idx              = inactive(idxMax);
        fl.max_inactive_flux = max_flux;
        fl.max_inactive_rxn  = rxn_idx;
        fprintf('  max inactive flux = %.3g on reaction %d\n', max_flux, rxn_idx);
    end

    % --- flags & validity as before ---
    fl.supp_size = numel(supp);
    fl.emptysupp = isempty(supp);
    if fl.emptysupp
        fl.nullity_ok = false;
        fl.balance_ok = false;
        fl.fluxpos_ok = false;
        fl.badpair_ok = false;
        fl.anchor_ok  = false;
        isValid = false;
        return;
    end

    tol_rank      = 1e-8;
    fl.nullity_ok = (numel(supp) - rank(full(S(:,supp)), tol_rank) == 1);
    fl.balance_ok = max(abs(S(:,supp) * v(supp))) <= tol_balance;
    fl.fluxpos_ok = all(abs(v(supp)) >= eps_flux);
    bp            = badPartner(supp);
    fl.badpair_ok = ~any(bp>0 & ismember(bp, supp));
    if isempty(anchorIdx)
        fl.anchor_ok = true;
    else
        fl.anchor_ok = ismember(anchorIdx, supp);
    end

    isValid = fl.nullity_ok && fl.balance_ok && fl.fluxpos_ok && fl.badpair_ok && fl.anchor_ok;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [isValid, supp] = oldcheckSupport(x, S, tol_balance, badPartner, eps_flux)
    % Extract flux and binary vectors
    n      = size(S,2);
    v      = x(1:n);
    b      = x(n+1:2*n);
    isValid = false;

    % 1) initial support from b
    supp = find(b > 0.5);
    if isempty(supp)
        return;
    end

    % 2) nullity=1 on original support
    tol_rank = 1e-8;
    if numel(supp) - rank(full(S(:,supp)), tol_rank) ~= 1
        disp('nullity');
        return;
    end

    % 3) flux lower-bound check
    if any(abs(v(supp)) < eps_flux)
        disp('flux too low');
        return;
    end

    % 4) steady-state residual on solver v
    residual = S(:, supp) * v(supp);
    res      = max(abs(residual));
    if res > tol_balance
        % 5) attempt greedy leak addition
        inactive = find(~ismember(1:n, supp));
        v_inactive = v(inactive);                   % solver fluxes on those reactions

        % compute the max absolute flux and its reaction index
        [ max_flux, idx ] = max(abs(v_inactive));
        rxn_idx = inactive(idx);

fprintf('  max inactive flux = %.3g on reaction %d\n', max_flux, rxn_idx);


        leaks    = inactive(abs(v(inactive)) > tol_balance);
        improved = true;

        while res > tol_balance && ~isempty(leaks) && improved
            improved = false;
            best_j   = [];
            best_res = res;

            % try each leak candidate
            for j = leaks
                cand  = [supp; j];
                % nullity check
                if numel(cand) - rank(full(S(:,cand)), tol_rank) ~= 1
                    disp("addition led to nullity not 1")
                    continue;
                end
                % residual on solver v
                res_j = max(abs(S(:,cand) * v(cand)));
                if res_j < best_res
                    best_j   = j;
                    best_res = res_j;
                end
            end

            % if one reaction improved residual, accept it
            if ~isempty(best_j)
                supp    = [supp; best_j];
                res     = best_res;
                leaks   = setdiff(leaks, best_j);
                improved = true;
                fprintf('  added leak %d, residual → %.2e\n', best_j, res);
            end
        end

        % if still out of tolerance => reject
        if res > tol_balance
            return;
        end
    end

    % 6) bad‐pair exclusion
    for i = supp'
        bp = badPartner(i);
        if bp>0 && any(supp==bp)
            return;
        end
    end

    % passed all checks
    isValid = true;
end
