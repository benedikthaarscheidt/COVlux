%% EFM_gurobi_pool_enumeration_optimized_nocard.m
% Fast enumeration of EFMs via Gurobi pool,
% with immediate no‐good cuts, vector bad‐pair test,
% seeding with retries, pooling without cardinality cuts.

%% 0) Gurobi threads & paths
nThreads = feature('numcores');
fprintf('Using %d threads for Gurobi\n',nThreads);
setenv('OMP_NUM_THREADS',num2str(nThreads));
addpath(fullfile(getenv('GUROBI_HOME'),'examples','matlab'),'-begin');

%% 1) Load irreversible pruned model
data     = load('/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/e_coli_core_splitandpruned.mat','pruned_ir');
model_ir = data.pruned_ir;
S        = model_ir.S;
rxnNames = model_ir.rxns;
[m,n]    = size(S);

%% 2) Parameters
eps_flux    = 9e-6;
tol_balance = 1e-6;     % steady‐state tolerance
M            = 1e3;
badPair_count = 0;
suppThreshZ       = 0.5;      % binary → support threshold
requireSuppBalance = true;    % also require Sv_supp * v_supp ≈ 0 (in addition to full Sv≈0)

%% 3) Build static MILP template
aeq = [S, sparse(m,n)];
beq = zeros(m,1);

A1 = [ speye(n), -M*speye(n)
       -speye(n),  eps_flux*speye(n) ];
b1 = zeros(2*n,1);

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
importE = find(importCols);
exportE = find(exportCols);

cutImp = sparse(1,2*n); cutImp(n+importE) = -1;
cutExp = sparse(1,2*n); cutExp(n+exportE) = -1;

Aineq = [A1; cutImp; cutExp];
bineq = [b1; -1; -1];

lb0    = [model_ir.lb; zeros(n,1)];
ub0    = [model_ir.ub; ones(n,1)];
vtype0 = [repmat('C',n,1); repmat('B',n,1)];

template.A          = [aeq; Aineq];
template.rhs        = [beq; bineq];
template.sense      = [ repmat('=',m,1); repmat('<',size(Aineq,1),1) ];
template.lb         = lb0;
template.ub         = ub0;
template.vtype      = vtype0;
template.modelsense = 'min';
template.obj        = [ zeros(n,1) ; ones(n,1) ];  % <<-- initialize objective

%% 4) Gurobi parameters
pPool = struct( ...
  'OutputFlag',    0, ...
  'PoolSearchMode',2, ...
  'PoolSolutions', 500, ...
  'TimeLimit',     600, ...
  'MIPGap',        0.0, ...
  'Threads',       nThreads );
pSingle = pPool;
pSingle.PoolSearchMode = 0;
pSingle.PoolSolutions  = 1;

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
maxSeedRetries = 40;
covered  = false(n,1);
skipped  = false(n,1);
allEFMs  = {};
allX     = {};
batch    = 0;
prevUc=0;



while any(~covered & ~skipped)
    batch = batch + 1;
    toTry = find(~covered & ~skipped);
    fprintf('\n=== Batch %d: %d reactions to try (%d/%d covered) ===\n', ...
            batch, numel(toTry), sum(covered), n);
    % pick next uncovered reaction
    
    % pick next uncovered reaction
    idx = randi(numel(toTry));
    i = toTry(idx);
    prevUc = i; 

    %% 6.1) Seeding with retries
    uc = toTry(1);
    seeded = false;
    for att = 1:maxSeedRetries
        M1          = template;           % local copy
        M1.lb(n+uc) = 1;                  % force z_uc = 1
        M1.obj      = template.obj;       % ensure obj is set

        sol1 = gurobi(M1, pPool);
        if ~ismember(sol1.status,{'OPTIMAL','TIME_LIMIT'})
            fprintf('  Seed infeasible for %d—retry %d/%d\n', uc,att,maxSeedRetries);
            continue;
        end

        x1    = getSolutionVector(sol1);
        v1    = x1(1:n);
        supp1=find(v1>=(eps_flux-tol_balance));
        %supp1 = find(x1(n+1:end)>0.5);

        % various validity checks
        nullOK    = (numel(supp1)-rank(S(:,supp1)))==1;
        fluxOK    = all(abs(v1(supp1))>=eps_flux);
        bpOK      = ~any(badPartner(supp1)>0 & ismember(badPartner(supp1),supp1));
        balanceOK = max(abs(S*v1)) <= tol_balance;    % <<-- steady‐state check
        if ~balanceOK
            disp("balance")
        end 
        if ~nullOK
            disp("nullity")
        end 
        if ~bpOK 
            [supp1, v_k] = cancel_twins(supp1, v1, badPartner, S, eps_flux, tol_balance);
        end 

        nullOK    = (numel(supp1)-rank(S(:,supp1)))==1;
        fluxOK    = all(abs(v1(supp1))>=eps_flux);
        bpOK      = ~any(badPartner(supp1)>0 & ismember(badPartner(supp1),supp1));
        balanceOK = max(abs(S*v1)) <= tol_balance;    % <<-- steady‐state check
        
        if ~balanceOK
            disp("balance")
        end 
        if ~nullOK
            disp("nullity")
        end 

        if nullOK && fluxOK && bpOK && balanceOK 
            fprintf('  Seed accepted (|supp|=%d)\n',numel(supp1));
            covered(supp1) = true;
            fprintf('    Coverage: %d/%d\n', sum(covered), n);
            allEFMs{end+1} = supp1;
            allX   {end+1} = v1;
            seeded = true;
            break;
        else
            fprintf('  Seed invalid—excluded (|supp|=%d)\n',numel(supp1));
        end
    end

    if ~seeded
        warning('  Skipping reaction %d after %d seeds\n', uc, maxSeedRetries);
        skipped(uc) = true;
        continue;
    end

    %% 6.2) Pooling (no cardinality)
    M2 = template;  % local copy for pooling

    % require at least one still‐uncovered reaction
    rem  = find(~covered);
    cutU = sparse(1,2*n); cutU(n+rem) = -1;
    M2.A     = [M2.A; cutU];
    M2.rhs   = [M2.rhs; -1];
    M2.sense = [M2.sense; '<'];
    M2.obj   = template.obj;             % <<-- ensure obj is set

    sol2 = gurobi(M2, pPool);
    poolSol = [];
    if isfield(sol2,'pool') && ~isempty(sol2.pool)
        poolSol = sol2.pool;
    elseif isfield(sol2,'xn')
        poolSol = struct('xn',sol2.xn);
    elseif isfield(sol2,'x')
        poolSol = struct('xn',sol2.x);
    end

    % post‐filter pooled candidates
    for p = 1:numel(poolSol)
        xk    = getSolutionVector(poolSol(p));
        v2    = xk(1:n);
        supp2 = find(xk(n+1:end)>0.6);

        if (numel(supp2)-rank(S(:,supp2)))~=1,                     continue; end
        if any(abs(v2(supp2))<eps_flux),                          continue; end
        if any(badPartner(supp2)>0 & ismember(badPartner(supp2),supp2)), continue; end
        if max(abs(S*v2))>tol_balance,                            continue; end  % steady‐state

        % accept and ban
        covered(supp2)   = true;
        allEFMs{end+1}   = supp2;
        allX   {end+1}   = v2;
        template         = applyNoGoodCut(template, supp2, n);
    end
end

%% 7) Report & save
stillUncov = find(~covered);
if ~isempty(stillUncov)
    fprintf('\nUncovered (%d): %s\n', numel(stillUncov), mat2str(stillUncov'));
else
    fprintf('\nAll covered.\n');
end
skList = find(skipped);
if ~isempty(skList)
    fprintf('Skipped (%d): %s\n', numel(skList), mat2str(skList'));
end
fprintf('\nDone. badPairs=%d\n', badPair_count);

% Build flux matrix & save
fluxMat = zeros(n,numel(allEFMs));
for i=1:numel(allEFMs)
    fluxMat(:,i) = allX{i};
end
save('EFMs_pool_badpair_nocard_3.mat','allEFMs','fluxMat','covered','skipped','badPair_count');
disp("fone")

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
    oth = setdiff(1:n,supp);
    cut = sparse(1,2*n);
    cut(n+supp) = 1;
    cut(n+oth)  = -1;
    T.A     = [T.A; cut];
    T.rhs   = [T.rhs; numel(supp)-1];
    T.sense = [T.sense; '<'];
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
                    %disp("the same")
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