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
%%

% 1) Decision variables x = [ v(1:n); b(1:n) ]
%    v_i continuous in [lb_i, ub_i], b_i binary
lb = [model_ir.lb; zeros(n,1)];
ub = [model_ir.ub; ones(n,1)];
ub_eff = model_ir.ub;
ub_eff(isinf(ub_eff)) = M;

vtype = [ repmat('C',n,1) ; repmat('B',n,1) ];

% 2) Mass balance:   S*v = 0
Aeq = [ sparse(S), sparse(m,n) ];
beq = zeros(m,1);

% 3) big‑M coupling & epsilon link:
%    v_i - M*b_i <= 0
%   -v_i + eps*b_i <= 0
A1 = [  speye(n),  -M*speye(n)
       -speye(n),  eps_flux*speye(n) ];
b1 = zeros(2*n,1);

% 4) import/export cuts:
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



% 5) bad‐pair bans: for each reversible pair (i,j),
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

% 6) assemble ALL inequalities
Aineq = [ A1
          Aimp
          Aexp ];
bineq = [ b1
          bimp
          bexp ];

sense = [ repmat('=', size(Aeq,1),1)
          repmat('<', size(Aineq,1),1) ];

% 7) objective: minimize sum(b) = [ zeros(1,n), ones(1,n) ] * x
obj = [ zeros(n,1) ; ones(n,1) ];

% 8) pack into Gurobi model struct
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
%technically we dont need the pooling settings here but I used them for
%playing around 
pPool = struct( ...
  'OutputFlag',     0, ...
  'PoolSearchMode',1, ...
  'PoolSolutions',100, ...
  'TimeLimit',    600, ...
  'MIPGap',0.5, ...
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
  'MIPGap',0.1, ...
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


fvaParams.OutputFlag   = 0;    % suppress solver output
fvaParams.TimeLimit    = 60;   % time limit per LP (seconds)
fvaParams.FeasibilityTol = 1e-9; % strict balance enforcement



%% 6) Enumerate 
covered = false(n,1);
allEFMs = cell(n,1);
allX    = cell(n,1);

maxRetries = 200;    

for i = 1:n
    if covered(i), continue; end
    fprintf("\n--- Reaction %d (%s) ---\n", i, rxnNames{i});

    localModel = linModel;        % start with full linearized MILP
    found = false;

    for attempt = 1:maxRetries
        % copy & force support reaction
        M1 = localModel;
        % force b_i = 1
        row = sparse(1,2*n); row(n+i) = 1;
        M1.A     = [M1.A; row];
        M1.rhs   = [M1.rhs; 1];
        M1.sense = [M1.sense; '='];

        % solve MILP to get a candidate support
        sol = gurobi(M1, pSeed);
        if ~isfield(sol,'status') || ~strcmp(sol.status,'OPTIMAL')
            fprintf('  attempt %d: no OPTIMAL (%s) → retrying\n', attempt, sol.status);
            continue;
        end
        x = sol.x;
        supp = find(x(n+1:2*n) > 0.5);

        % 1) nullity check
        [null_dim, ~, ~] = nullity_and_tol(S(:, supp));
        if null_dim ~= 1
            % truly non-EFM: ban permanently
            localModel = applyNoGoodCut(localModel, supp, n);
            continue;
        end

        % 2) enforce steady-state by LP on same support
        % build small LP: S(:,supp)*v_supp = 0, eps <= v_supp <= ub_eff(supp)
        LP.A        = sparse(S(:, supp));
        LP.rhs      = zeros(m,1);
        LP.sense    = repmat('=', m,1);
        LP.lb       = eps_flux * ones(numel(supp),1);
        LP.ub       = ub_eff(supp);
        LP.vtype    = repmat('C', numel(supp),1);
        LP.modelsense = 'min';
        LP.obj      = zeros(numel(supp),1);

        solLP = gurobi(LP, fvaParams);
        if isfield(solLP,'status') && strcmp(solLP.status,'OPTIMAL')
            % steady-state feasible: accept
            v_full = zeros(n,1);
            v_full(supp) = solLP.x;
            allEFMs{i} = supp;
            allX{i}    = v_full;
            covered(supp) = true;
            fprintf('   → valid EFM (|supp|=%d) on attempt %d\n', numel(supp), attempt);
            found = true;
            break;
        else
            % false alarm: support cannot meet steady state -> ban
            localModel = applyNoGoodCut(localModel, supp, n);
            continue;
        end
    end

    if ~found
        warning('   No valid EFM for reaction %d after %d tries', i, maxRetries);
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

function T = applyNoGoodCut(T, supp, n)
    oth = setdiff(1:n, supp);
    cut = sparse(1,2*n);
    cut(n+supp) = 1;
    cut(n+oth)  = -1;
    T.A     = [T.A; cut];
    T.rhs   = [T.rhs; numel(supp)-1];
    T.sense = [T.sense; '<'];
end


function [isValid, supp] = checkSupport(x, S, tol_balance, badPartner, eps_flux)
    n    = size(S,2);
    v      = x(1:n);
    b      = x(n+1:2*n);
    isValid = false;

    % 1) initial support from b
    supp = find(b > 0.5);
    %supp    = find(abs(v) >= eps_flux);
    if isempty(supp), return; end
    tol_rank = 1e-8;
    if numel(supp) - rank(S(:,supp), tol_rank) ~= 1
        disp("nullity")
        return;
    end
    if any(abs(v(supp)) < (eps_flux-tol_balance))
        disp("flux")
        return;
    end
    if max(abs(S(:,supp)*v(supp))) > tol_balance
        disp("steady")
        return;
    end
    if any(badPartner(supp)>0 & ismember(badPartner(supp), supp))
        return;
    end
    isValid = true;
end


default_tol_rank = 1e-8;
function [null_dim, singvals, tol_rank] = nullity_and_tol(mat)
    % full SVD
    singvals = svd(mat, 'econ');
    max_sv = singvals(1);
    tol_rank = max(size(mat)) * max_sv * eps;
    rank_est = sum(singvals > tol_rank);
    null_dim = size(mat,2) - rank_est;
end
