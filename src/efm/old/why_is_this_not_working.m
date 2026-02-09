%% 0) Gurobi threads & paths
nThreads = feature('numcores');
fprintf('Using %d threads for Gurobi\n',nThreads);
setenv('OMP_NUM_THREADS',num2str(nThreads));
addpath( fullfile(getenv('GUROBI_HOME'),'examples','matlab'), '-begin' );

%% 1) Load irreversible pruned model
data     = load('/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/e_coli_core_splitandpruned.mat','pruned_ir');%e_coli_core_splitandpruned%/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/pruned_iML1515_irrev_split_biomassrm_v5.mat'
model_ir = data.pruned_ir;
S        = model_ir.S;
rxnNames = model_ir.rxns;
[m,n]    = size(S);
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
%% 2) Parameters
eps_flux      = 1e-5;
tol_bal       = 1e-7;
M             = 1e3;
badPair_count = 0;

%% 3) Precompute bad‐pair lookup
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

%% 4) constraint matrix

%note that with the badpair constraint the enumeration is not possible for
%some reactions. I have commented out these constraints to run a less
%constraint model and exclude bad pairs later. 
iv = 1:n;
ib = n+(1:n);
iy = 2*n+(1:n);

Abad = sparse(badPair_count, 3*n);
bbad = ones(badPair_count, 1);
k = 1;
for i = 1:n
    j = badPartner(i);
    if j > 0 && i < j
        Abad(k, n+i) = 1;
        Abad(k, n+j) = 1;
        k = k + 1;
    end
end

aeq1 = [S,            sparse(m,n), sparse(m,n)];
beq1 = zeros(m,1);
aeq2 = [sparse(n,n),  speye(n),    speye(n)];
beq2 = ones(n,1);
aeq  = [aeq1; aeq2];
beq  = [beq1; beq2];

A1  = [-speye(n),  eps_flux*speye(n),  sparse(n,n)];
b1  = zeros(n,1);

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
cutImp = sparse(1,3*n); cutImp(n+find(importCols)) = -1;
cutExp = sparse(1,3*n); cutExp(n+find(exportCols)) = -1;

Aineq = [A1; cutImp; cutExp]; % append Abad here if desired
bineq = [b1; -1; -1];         % append bbad if Abad is used

lb0    = [model_ir.lb; zeros(n,1); zeros(n,1)];
ub0    = [model_ir.ub; ones(n,1);  ones(n,1)];
vtype0 = [repmat('C',n,1); repmat('B',n,1); repmat('C',n,1)];

bareTemplate.A          = [aeq; Aineq];
bareTemplate.rhs        = [beq; bineq];
bareTemplate.sense      = [repmat('=',size(aeq,1),1); repmat('<',size(Aineq,1),1)];
bareTemplate.lb         = lb0;
bareTemplate.ub         = ub0;
bareTemplate.vtype      = vtype0;
bareTemplate.modelsense = 'min';
bareTemplate.obj        = [zeros(n,1); ones(n,1); zeros(n,1)];

for j = 1:n
    bareTemplate.sos(j).type   = 1;
    bareTemplate.sos(j).index  = [iv(j), iy(j)];
    bareTemplate.sos(j).weight = [1, 2];
end


%% 4) Gurobi parameters
commonParams = struct( ...
   'FeasibilityTol',1e-9, ...
   'IntFeasTol',    1e-9, ...
   'NumericFocus',  3);
%technically we dont need the pooling settings here but I used them for
%playing around 
pPool = struct( ...
  'OutputFlag',     0, ...
  'PoolSearchMode',1, ...
  'PoolSolutions',100, ...
  'TimeLimit',    600, ...
  'MIPGap',0.0, ...
  'Cuts',0, ...
  'Presolve',0, ...%this was 2 
  'Heuristics',0.5, ...
  'MIPFocus',1, ...
  'Threads',nThreads );

pSeed = struct( ...
  'OutputFlag',     0, ...
  'PoolSearchMode',0, ... %one optimal solution
  'PoolSolutions',1, ...
  'TimeLimit',    600, ...
  'MIPGap',0.0, ...
  'Cuts',-1, ... %automatic cut determining 
  'Presolve',0, ...
  'Heuristics',0.1, ...
  'MIPFocus',3, ... %aggressive MILP solve 
  'Threads',nThreads );
for p = {pPool,pSeed}
    p{1}.FeasibilityTol = commonParams.FeasibilityTol;
    p{1}.IntFeasTol     = commonParams.IntFeasTol;
    p{1}.NumericFocus   = commonParams.NumericFocus;
end




%% 6) Enumerate
covered = false(n,1);
allEFMs = cell(n,1);
allX    = cell(n,1);
maxRetries = 200;

for i = 1:n
    if covered(i), continue; end
    fprintf('\n--- Reaction %d (%s) ---\n', i, rxnNames{i});
    localT = bareTemplate;
    found = false;
    for attempt = 1:maxRetries
        M1 = localT;
        cutF = sparse(1,3*n); cutF(n+i) = 1;
        M1.A     = [M1.A;   cutF];
        M1.rhs   = [M1.rhs; 1];
        M1.sense = [M1.sense; '>'];
        M1.lb(n+i) = 1;
        M1.ub(n+i) = 1;
        
        sanityCheckSOSIndices(M1, n);

        disp(M1.sos(1));
        sol = gurobi(M1, pSeed);
        if ~isfield(sol,'status') || ~strcmp(sol.status,'OPTIMAL')
            fprintf('  attempt %d: no OPTIMAL (%s)\n', attempt, sol.status);
            continue;
        end

        x = getSolutionVectorMILP(sol);
        assertLinkingSOS(x, eps_flux, 1e-9);
        v = x(1:n);

        [ok, supp] = checkSupport_bounds(x, S, model_ir.lb, model_ir.ub, tol_bal, eps_flux, badPartner);
        if ok
            fprintf('   → valid EFM (|supp|=%d) on attempt %d\n', numel(supp), attempt);
            covered(supp) = true;
            allEFMs{i} = supp;
            allX{i}    = v;
            found = true;
            break;
        else
            fprintf('   attempt %d: invalid support (|supp|=%d)\n', attempt, numel(supp));
            localT = applyNoGoodCut3n(localT, supp, n);
        end
    end
    if ~found
        warning('   Could not find a valid EFM for reaction %d after %d tries\n', i, maxRetries);
    end
end
%% 7) Collect & save 
EFMs_found = allEFMs(~cellfun('isempty',allEFMs));
fluxMat    = zeros(n, numel(EFMs_found));
for k = 1:numel(EFMs_found)
    fluxMat(:,k) = allX{ find(~cellfun('isempty',allEFMs),1,'first') + (k-1) };
end
save('EFMs_retry_per_rxn.mat','EFMs_found','fluxMat','badPair_count');
fprintf('Done: recovered %d EFMs (with up to %d retries each)\n', numel(EFMs_found), maxRetries);


%% --- helper functions ---
function x = getSolutionVector(sol)
% Pull out the primal x from sol.x or sol.xn or sol.pool(1).xn
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

function x = getSolutionVectorMILP(sol)
    if isfield(sol,'x') && ~isempty(sol.x)
        x = sol.x;
    else
        error('MILP incumbent (sol.x) missing; refusing to use sol.xn/pool.');
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


function [isValid, supp] = checkSupport_bounds(x, S, lb, ub, tol_balance, eps_flux, badPartner)
    n = size(S,2);
    v = x(1:n);
    b = x(n+1:2*n);
    isValid = false;
        thr = max(0, eps_flux - tol_balance);
    supp = find(v >= thr);
    
    if any(b(supp) < 0.5) || any(v(b>0.5) < thr)
        warning('b and v are inconsistent → tighten M/eps or clean support by flux.');
    end
    if isempty(supp), return; end
    r = S(:,supp)*v(supp);
    if norm(r,Inf) > tol_balance
        disp("steady")
        return;
    end
    A = full(S(:,supp));
    if ~isempty(ub)
        actUB = find(abs(v(supp) - ub(supp)) <= 1e-9);
        if ~isempty(actUB)
            E = sparse(1:numel(actUB), actUB, 1, numel(actUB), numel(supp));
            A = [A; full(E)];
        end
    end
    s = svd(A,'econ');
    thr = max(size(A))*eps(max(s));
    rk = sum(s > max(thr,1e-10));
    if numel(supp) - rk ~= 1
        disp("nullity") 
        return; 
    end
    if any(badPartner(supp)>0 & ismember(badPartner(supp), supp)), return; end
    isValid = true;
end

function [ok, w_full] = refit_on_support_nonneg(S, lb, ub, supp, eps_flux)
    m = size(S,1); n = size(S,2);
    w_full = zeros(n,1);
    if isempty(supp), ok=false; return; end
    Aeq = [S(:,supp), sparse(m,1); ones(1,numel(supp)), 0];
    beq = [zeros(m,1); 1];
    lbA = max(lb(supp), 0) + eps_flux;
    ubA = ub(supp);
    model.A = sparse(Aeq);
    model.rhs = beq;
    model.sense = [repmat('=',m+1,1)];
    model.lb = [lbA; 1];
    model.ub = [ubA; 1];
    model.vtype = repmat('C',numel(supp)+1,1);
    model.modelsense = 'min';
    model.obj = zeros(numel(supp)+1,1);
    sol = gurobi(model, struct('OutputFlag',0,'FeasibilityTol',1e-9,'IntFeasTol',1e-9));
    if isfield(sol,'status') && strcmp(sol.status,'OPTIMAL')
        w = sol.x(1:numel(supp));
        w_full(supp) = w;
        ok = true;
    else
        ok = false;
    end
end


function assertLinking(x, eps_flux, M, feasTol)
    n = numel(x)/2;
    v = x(1:n); b = x(n+1:2*n);
    assert(all(v - M*b <= feasTol + 1e-12), 'Linking violated: v - M b <= 0');
    assert(all(-v + eps_flux*b <= feasTol + 1e-12), 'Linking violated: -v + eps b <= 0');
end



function T = applyNoGoodCut3n(T, supp, n)
    oth = setdiff(1:n, supp);
    cut = sparse(1,3*n);
    cut(n+supp) = 1;
    cut(n+oth)  = -1;
    T.A     = [T.A; cut];
    T.rhs   = [T.rhs; numel(supp)-1];
    T.sense = [T.sense; '<'];
end

function assertLinkingSOS(x, eps_flux, feasTol)
    n = numel(x)/3;
    v = x(1:n);
    b = x(n+1:2*n);
    y = x(2*n+1:3*n);
    
    tolEq  = max(10*feasTol, 1e-8);   % equality/coupling tolerance
    tolSOS = max(50*feasTol, 1e-7);   % SOS overlap tolerance
    
    if max(abs((b + y) - 1)) > tolEq
        error('b + y = 1 violated (max viol %.3e).', max(abs((b+y)-1)));
    end
    if any(-v + eps_flux*b > tolEq)
        error('-v + eps*b <= 0 violated (max viol %.3e).', max(-v + eps_flux*b));
    end
    
    mask = (v > tolSOS) & (y > tolSOS);
    if any(mask)
        j = find(mask,1,'first');
        error('SOS1 {v_j,y_j} violated at j=%d: v=%.3e, y=%.3e, b=%.3f', j, v(j), y(j), b(j));
    end
    
    if any(v < -tolEq)
        error('v must be nonnegative (min v = %.3e).', min(v));
    end
end

function sanityCheckSOSIndices(M, n)
  assert(isfield(M,'sos') && numel(M.sos)==n, 'SOS blocks missing/wrong length.');
  for j=1:n
    idx = M.sos(j).index;
    if numel(idx)~=2 || idx(1)~=j || idx(2)~=2*n+j
      error('SOS index mismatch at j=%d: [%d %d] vs [%d %d].', j, idx(1), idx(2), j, 2*n+j);
    end
    if M.sos(j).type~=1, error('SOS at j=%d not type 1.', j); end
  end
end 