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
Abad = sparse(badPair_count, 2*n);
bbad = ones(badPair_count, 1);
k = 1;
for i = 1:n
    j = badPartner(i);
    if j > 0 && i < j  % b_i + b_j <= 1 for each pair
        Abad(k, n+i) = 1;
        Abad(k, n+j) = 1;
        k = k + 1;
    end
end


aeq = [S, sparse(m,n)];
beq = zeros(m,1);
A1  = [ speye(n), -M*speye(n)
       -speye(n),  eps_flux*speye(n) ];
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

% Add bad pair constraints to the inequality system
Aineq = [A1; cutImp; cutExp];% Abad];
bineq = [b1; -1; -1];% bbad];

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




%% 6) Enumerate 
covered = false(n,1);
allEFMs = cell(n,1);
allX    = cell(n,1);

maxRetries = 200;    

for i = 1:n
    if covered(i)
        continue;
    end
    fprintf('\n--- Reaction %d (%s) ---\n', i, rxnNames{i});

    localT = bareTemplate;
    
    found = false;
    for attempt = 1:maxRetries
        % 1) force reaction i active
        
        M1 = localT;
        cutF = sparse(1,2*n); cutF(n+i) = 1;
        M1.A      = [M1.A;   cutF];
        M1.rhs    = [M1.rhs; 1];
        M1.sense  = [M1.sense; '>'];
        M1.lb(n+i) = 1;
        M1.ub(n+i) = 1;
        %flux_val = range(i,1) + (range(i,2) - range(i,1)) * rand(1,1);
        %M1.lb(i)   = flux_val; %Zoran suggested that this is going to solve the problem with the steady state but it does not
        %M1.ub(i)   = flux_val; %Zoran suggested that this is going to solve the problem with the steady state but it does not
       

        % 2) solve
        sol = gurobi(M1, pSeed);
        
        if ~isfield(sol,'status') || ~strcmp(sol.status,'OPTIMAL')
            fprintf('  attempt %d: no OPTIMAL (%s) → retrying\n', attempt, sol.status);
            continue;
        end
        % 3) extract and test support
        x    = getSolutionVector(sol);
        v    = x(1:n);
        z    = x(n+1:end);
        

        [ok, supp] = checkSupport(x, S, tol_bal, badPartner, eps_flux);

        if ok
            % valid EFM found
            fprintf('   → valid EFM (|supp|=%d) on attempt %d\n', numel(supp), attempt);
            covered(supp)   = true;
            allEFMs{i}      = supp;
            allX   {i}      = v;
            found = true;
            break;
        else
            %fprintf('   attempt %d: invalid support (|supp|=%d) → banning and retrying\n', attempt, numel(supp));
            localT = applyNoGoodCut(localT, supp, n); %we only need to do it for the local template as when we break this loop a new, uncovered reaction will be chosen which guarantees a different support.
        end
    end

    if ~found
        warning('   Could not find a valid EFM for reaction %d after %d tries\n', i, maxRetries);
    end
    % — note: we do *not* apply a no‑good cut for the valid one, since
    %   forcing z_{i+1}=1 next iteration guarantees a different support.
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
    if max(abs(S*v)) > tol_balance
        disp("steady")
        return;
    end
    inactive = find(~ismember(1:n, supp));
    v_inactive = v(inactive);                   % solver fluxes on those reactions

    % compute the max absolute flux and its reaction index
    [max_flux, idx ] = max(abs(v_inactive));
    rxn_idx = inactive(idx);
    fprintf('  max inactive flux = %.3g on reaction %d\n', max_flux, rxn_idx);
    if any(badPartner(supp)>0 & ismember(badPartner(supp), supp))
        return;
    end
    isValid = true;
end
