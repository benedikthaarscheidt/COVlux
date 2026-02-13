% with immediate no‐good cuts for invalid supports only,
% vector bad‐pair test, seeding with retries, pooling without cardinality cuts

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
eps_flux      = 1e-5;
tol_balance   = 1e-6;
M             = 1e3;
badPair_count = 0;

%% 3) Build static MILP template with indicator constraints
aeq = [S, sparse(m,n)];
beq = zeros(m,1);
ub_eff = model_ir.ub;
ub_eff(isinf(ub_eff)) = M;

% Initialize inequality constraints (will add indicator constraints separately)
Aineq = [];
bineq = [];

% import/export cuts (same as before)
nz = (S~=0);
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

Aineq = [Aineq; cutImp; cutExp];
bineq = [bineq; -1; -1];

% Variable bounds and types (same as before)
lb0 = [model_ir.lb; zeros(n,1)];
ub0 = [ub_eff; ones(n,1)];
vtype0 = [repmat('C',n,1); repmat('B',n,1)];

% Create model structure
bareTemplate = struct();
bareTemplate.A = [aeq; Aineq];
bareTemplate.rhs = [beq; bineq];
bareTemplate.sense = [repmat('=',m,1); repmat('<',size(Aineq,1),1)];
bareTemplate.lb = lb0;
bareTemplate.ub = ub0;
bareTemplate.vtype = vtype0;
bareTemplate.modelsense = 'min';
bareTemplate.obj = [zeros(n,1); ones(n,1)];


%% Add indicator constraints (replaces your A1/b1 constraints)
bareTemplate.genconind = repmat(struct(...
    'binvar',0, 'binval',0, 'a',[], 'sense','', 'rhs',0), 2*n, 1);

for j = 1:n
    % b_j = 0 ⇒ v_j = 0
    bareTemplate.genconind(2*j-1).binvar = n + j;
    bareTemplate.genconind(2*j-1).binval = 0;
    bareTemplate.genconind(2*j-1).a = sparse(j,1,1,2*n,1);
    bareTemplate.genconind(2*j-1).sense = '=';
    bareTemplate.genconind(2*j-1).rhs = 0;
    
    % b_j = 1 ⇒ v_j ≥ eps_flux
    bareTemplate.genconind(2*j).binvar = n + j;
    bareTemplate.genconind(2*j).binval = 1;
    bareTemplate.genconind(2*j).a = sparse(j,1,1,2*n,1);
    bareTemplate.genconind(2*j).sense = '>';
    bareTemplate.genconind(2*j).rhs = eps_flux;
end

%% 4) Gurobi parameters
commonParams = struct( ...
    'FeasibilityTol', 1e-9, ...
    'IntFeasTol', 1e-9, ...
    'NumericFocus', 3, ... % Highest numerical stability
    'ScaleFlag', 2, ...   % Aggressive scaling
    'Quad', 1);           % For numerical issues
pPool = struct( ...
  'OutputFlag',     0, ...
  'PoolSearchMode',2, ...
  'PoolSolutions',100, ...
  'TimeLimit',    60, ...
  'MIPGap',0.5, ...
  'Cuts',1, ...
  'Presolve',2, ...%this was 2 
  'Heuristics',0.5, ...
  'MIPFocus',1, ...
  'Threads',nThreads );
pSeed = struct( ... 
  'OutputFlag',     0, ...
  'PoolSearchMode',1, ...
  'PoolSolutions',1, ...
  'TimeLimit',    600, ...
  'MIPGap',0.5, ...
  'Cuts',1, ... 
  'Presolve',2, ...
  'Heuristics',0.5, ...
  'MIPFocus',2, ...
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


%% 


maxSeedAttempts = 50;
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
            %for j=[1:10]   % sample a few reactions
            %    fprintf('rxn %d: b=%g, v=%g\n', j, xk(n+j), xk(j));
            %end
            [ok, supp] = checkSupport(xk, S, tol_balance, badPartner, eps_flux);
            suppList{k} = supp;
            if ok && ismember(uc, supp)
                valid(k) = true;
                %vList{k} = xk(1:n);
                vList{k} = xk(supp);
            end
        end

        % ban only invalid supports on localTemplate
        for k = find(~valid).'
            localTemplate = applyNoGoodCut(localTemplate, suppList{k}, n);
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

    % POOLING: restart from globalTemplate (only permanent cuts)
    %localTemplate = globalTemplate;
    M2 = localTemplate;
    rem = find(~covered);
    cutU = sparse(1,2*n); cutU(n+rem) = 1;
    M2.A      = [M2.A;   cutU];
    M2.rhs    = [M2.rhs; 1];
    M2.sense  = [M2.sense; '>'];
    M2.lb(n+uc) = 1;
    sol2 = gurobi(M2, pPool);

     if isfield(sol2,'pool') && ~isempty(sol2.pool)
        poolSol = sol2.pool;
    elseif isfield(sol2,'xn')
        poolSol = struct('xn',sol2.xn);
    else
        poolSol = struct('xn',sol2.x);
    end

    numPools = numel(poolSol);
    valid    = false(numPools,1);
    suppList = cell(numPools,1);
    vList    = cell(numPools,1);
    parfor p = 1:numPools
        xk = getSolutionVector(poolSol(p));
        [ok, supp] = checkSupport(xk, S, tol_balance, badPartner, eps_flux);
        suppList{p} = supp;
        if ok
            valid(p) = true;
            %vList{p} = xk(1:n);
            vList{p} = xk(supp);
        end
    end

    % ban invalid pool supports and accept valid ones
    for p = 1:numPools
        if valid(p)
            covered(suppList{p})      = true;
            allEFMs{end+1}            = suppList{p};
            allX{end+1}               = vList{p};
            %globalTemplate            = applyNoGoodCut(globalTemplate, suppList{p}, n);
        end
    end
    fprintf('  Pooling: accepted %d modes\n', sum(valid));
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
fprintf('\nDone. badPairs=%d\n', badPair_count);

fluxMat = zeros(n,numel(allEFMs));
for i = 1:numel(allEFMs)
    fluxMat(:,i) = allX{i};
end
save('EFMs_pool_badpair_nocard_parpool_2.mat','allEFMs','fluxMat','covered','skipped','badPair_count');
disp("fone");

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

function [isValid, supp] = checkSupport(x, S, tol_balance, badPartner, eps_flux)
    isValid=false;
    n       = size(S,2);
    v       = x(1:n);
    b       = x(n+1:2*n);
    
    % 1) support from exact binaries
    supp    = find(b > 0.999);
    if isempty(supp)
        return;
    end
    inactive = setdiff(1:n, supp);
    maxInactive = max(abs(v(inactive)));
    minActive   = min(v(supp));
    nullity=numel(supp) - rank(S(:,supp), 1e-8) ;
    fprintf('  maxInactive=%g  minActive=%g  nullity=%d  res=%g\n',       maxInactive, minActive, nullity, max(abs(S*v)));
    % 2) correct nullity = 1
    if nullity ~= 1
        disp("nullity not 1 ")
        return;
    end

    % 3) flux‐value checks
    
    if any( v(supp) < eps_flux  - tol_balance)  
        disp("flux too small")
        disp(v(supp))% active too small
        return;
    end
    %if any( abs(v(inactive)) > tol_balance )         % inactive not zero
    %    return;
    %end

    % 4) mass‐balance residual
    if max( abs( S(:,supp) * v(supp) ) ) > tol_balance
        disp("steady state")
        return;
    end

    % 5) bad‐pair exclusion
    if any( badPartner(supp)>0 & ismember(badPartner(supp), supp) )
        return;
    end

    % if we made it here, it’s a valid EFM support
    isValid = true;
end
