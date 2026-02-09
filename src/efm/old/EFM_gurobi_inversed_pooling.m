%% EFM_gurobi_pool_enumeration_full_with_badpair.m
% Fast enumeration of EFMs via Gurobi pool,
% with dynamic split‐exclusivity, import/export, min‐flux, seeding with retries,
% pooling, bad‐pair counting, and no‐good cuts.

%% 0) Gurobi threads & paths
nThreads = feature('numcores');
fprintf('Detected %d CPU threads; using multi‐threaded Gurobi\n', nThreads);
setenv('OMP_NUM_THREADS', num2str(nThreads));
gurobiFolder = fullfile(getenv('GUROBI_HOME'),'examples','matlab');
addpath(gurobiFolder,'-begin');

%% 1) Load irreversible pruned model
data     = load('/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/e_coli_core_splitandpruned.mat','pruned_ir');%load('/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/pruned_iML1515_irrev_split_biomassrm_v4.mat','pruned_ir');
model_ir = data.pruned_ir;
[m,n]    = size(model_ir.S);
S        = model_ir.S;

%%


%% 2) Parameters
eps_flux      = 1e-6;
M             = 1e4;
badPair_count = 0;

%% 3) Build static MILP (template)
% 3.1) Mass‐balance
aeq = [model_ir.S, sparse(m,n)];
beq = zeros(m,1);

% 3.2) v–z coupling: v ≤ M z, v ≥ ε z
Aineq = [ ...
    speye(n),     -M*speye(n);      % v_i ≤ M z_i
   -speye(n),      eps_flux*speye(n) % v_i ≥ ε z_i
];
bineq = [ zeros(n,1)
          zeros(n,1) ];

% 3.3) Import/export requirement
exCols  = find(arrayfun(@(j) any(model_ir.S(:,j)~=0) && ...
    (all(model_ir.S(model_ir.S(:,j)~=0,j)>0) || ...
     all(model_ir.S(model_ir.S(:,j)~=0,j)<0)), 1:n));
importE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)>0), exCols));
exportE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)<0), exCols));

cutImp = sparse(1,2*n);  cutImp(   n+importE) = -1;  % -∑ z(import) ≤ -1
cutExp = sparse(1,2*n);  cutExp(   n+exportE) = -1;  % -∑ z(export) ≤ -1

Aineq = [Aineq; cutImp; cutExp];
bineq = [bineq; -1; -1];

% 3.4) Bounds & integer types
lb0    = [model_ir.lb; zeros(n,1)];
ub0    = [model_ir.ub; ones(n,1)];
vtype0 = [repmat('C',n,1); repmat('B',n,1)];

% 3.5) Assemble Gurobi template
template.A          = [aeq; Aineq];
template.rhs        = [beq; bineq];
template.sense      = [repmat('=',m,1); repmat('<',size(Aineq,1),1)];
template.lb         = lb0;
template.ub         = ub0;
template.vtype      = vtype0;
template.modelsense = 'min';

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
  'MIPGap',0.5, ...
  'Cuts',1, ...
  'Presolve',2, ...%this was 2 
  'Heuristics',0.5, ...
  'MIPFocus',1, ...
  'Threads',nThreads );
pSeed = struct( ...
  'OutputFlag',     0, ...
  'PoolSearchMode',1, ...
  'PoolSolutions',40, ...
  'TimeLimit',    600, ...
  'MIPGap',0.4, ...
  'Cuts',1, ... 
  'Presolve',2, ...
  'Heuristics',0.2, ...
  'MIPFocus',2, ...
  'Threads',nThreads );
for p = {pPool,pSeed}
    p{1}.FeasibilityTol = commonParams.FeasibilityTol;
    p{1}.IntFeasTol     = commonParams.IntFeasTol;
    p{1}.NumericFocus   = commonParams.NumericFocus;
end


%% 5) Enumeration loop with seed retry and final reporting
maxSeedRetries = 30;          % max seeding attempts per reaction
covered  = false(n,1);       % reactions covered by accepted EFMs
skipped  = false(n,1);       % reactions skipped after retries
allEFMs  = {};               % cell array of EFM supports
batch    = 0;

maxSeedRetries = 30;                % max seeding attempts per reaction
covered  = false(n,1);             % reactions covered by accepted EFMs
skipped  = false(n,1);             % reactions skipped after retries
allEFMs  = {};                     % cell array of EFM supports
allX     = {};                     % cell array of corresponding flux vectors
batch    = 0;

while any(~covered & ~skipped)
    batch = batch + 1;
    toTry = find(~covered & ~skipped);
    fprintf("\n=== Batch %d: %d reactions to try ===\n", batch, numel(toTry));

    %% 5.1) Phase I: seeding with retry limit
    uc = toTry(1);
    seedRetries = 0;
    seeded      = false;
    while seedRetries < maxSeedRetries
        % build & solve MILP forcing reaction uc on
        M1          = template;
        M1.lb(n+uc) = 1;
        M1.obj      = [zeros(n,1); ones(n,1)];
        sol1        = gurobi(M1, pSeed);

        if ismember(sol1.status, {'OPTIMAL','TIME_LIMIT'})
            x1    = getSolutionVector(sol1);
            supp1 = find(x1(n+1:end) > 0.5);

            % checks: nullity, min-flux, bad-pair
            if numel(supp1)-rank(model_ir.S(:,supp1))==1 && ...
               all(abs(x1(supp1)) >= eps_flux)          && ...
               ~hasBadPair(supp1, model_ir.rxns) && all(abs(S(:,supp1) * x1(supp1)) < 1e-7)

                fprintf("  Seed found EFM for reaction %d (|supp|=%d)\n", uc, numel(supp1));
                % record support and flux
                allEFMs{end+1} = supp1;
                allX   {end+1} = x1(1:n);
                covered(supp1) = true;
                template       = applyNoGoodCut(template, supp1, n);
                seeded         = true;
                break;
            else
                fprintf("  Seed candidate invalid—adding no-good cut and retrying\n");
                template = applyNoGoodCut(template, supp1, n);
            end
        else
            fprintf("  Seed infeasible for %d—retry %d/%d\n", uc, seedRetries+1, maxSeedRetries);
        end
        seedRetries = seedRetries + 1;
    end

    if ~seeded
        warning("  Could not cover reaction %d after %d seeds—skipping.\n", uc, maxSeedRetries);
        skipped(uc) = true;
        continue;
    end

    %% 5.2) Phase II: pooling additional EFMs of same size
    k  = numel(supp1);
    M2 = template;
    % enforce exact cardinality: sum(z) == k
    cutK_ub = sparse(1,2*n); cutK_ub(n+1:2*n) =  1;
    cutK_lb = sparse(1,2*n); cutK_lb(n+1:2*n) = -1;
    M2.A     = [M2.A; cutK_ub; cutK_lb];
    M2.rhs   = [M2.rhs; k; -k];
    M2.sense = [M2.sense; '<'; '>'];
    % require at least one still-uncovered reaction
    uc2 = find(~covered);
    cutU = sparse(1,2*n); cutU(n+uc2) = -1;
    M2.A     = [M2.A; cutU];
    M2.rhs   = [M2.rhs; -1];
    M2.sense = [M2.sense; '<'];
    M2.obj   = [zeros(n,1); ones(n,1)];

    sol2 = gurobi(M2, pPool);
    % extract pool solutions robustly
    if isfield(sol2,'pool') && ~isempty(sol2.pool)
        poolSol = sol2.pool;
    elseif isfield(sol2,'xn')
        poolSol = struct('xn', sol2.xn);
    elseif isfield(sol2,'x')
        poolSol = struct('xn', sol2.x);
    else
        warning('  No solution vector found in this pool iteration—skipping pooling phase.');
        poolSol = [];
    end

    % post-filter each pooled candidate
    for p = 1:numel(poolSol)
        xk   = getSolutionVector(poolSol(p));
        supp = find(xk(n+1:end) > 0.5);
        if numel(supp)-rank(model_ir.S(:,supp))~=1,        continue; end
        if any(abs(xk(supp)) < eps_flux),                  continue; end
        if hasBadPair(supp, model_ir.rxns),                continue; end
        if any(abs(S(:,supp1) * x1(supp1)) > 1e-7)         continue; end

        fprintf("  Pool accepted EFM (|supp|=%d)\n", numel(supp));
        allEFMs{end+1} = supp;
        allX   {end+1} = xk(1:n);
        covered(supp)  = true;
        template       = applyNoGoodCut(template, supp, n);
    end
end

%% % 6.1) Report coverage & skips
stillUncovered = find(~covered);
if ~isempty(stillUncovered)
    fprintf('\nUncovered reactions (%d):\n', numel(stillUncovered));
    for idx = stillUncovered(:)'
        fprintf('  %4d: %s\n', idx, model_ir.rxns{idx});
    end
else
    fprintf('\nAll reactions covered.\n');
end

skippedList = find(skipped);
if ~isempty(skippedList)
    fprintf('\nReactions skipped after %d failed seeds:\n', maxSeedRetries);
    for idx = skippedList(:)'
        fprintf('  %4d: %s\n', idx, model_ir.rxns{idx});
    end
end

fprintf('\nDone… badPairs=%d\n', badPair_count);

% 6.2) Build and save the flux matrix for the initial enumeration
% (Assumes you collected allX during enumeration)
numEFMs_original = numel(allEFMs);
fluxMat_original = zeros(n, numEFMs_original);
for k = 1:numEFMs_original
    fluxMat_original(:,k) = allX{k};
end

% Save initial results
save('EFMs_pool_badpair.mat', ...
     'allEFMs', 'fluxMat_original', 'covered', 'badPair_count');

% Also write CSV of the initial flux matrix
T0 = array2table(fluxMat_original, ...
    'RowNames',      model_ir.rxns, ...
    'VariableNames', strcat("EFM", string(1:numEFMs_original)) );
writetable(T0, 'EFMs_flux_initial.csv', 'WriteRowNames', true);

%% 6.3) Load and rescue
data = load('EFMs_pool_badpair.mat', 'allEFMs', 'fluxMat_original', 'covered', 'badPair_count');
allEFMs_original     = data.allEFMs;
fluxMat_original     = data.fluxMat_original;
covered_existing     = data.covered;
badPair_count        = data.badPair_count;


% Recompute boundaryE
exCols  = find(arrayfun(@(j) any(model_ir.S(:,j)~=0) && ...
    (all(model_ir.S(model_ir.S(:,j)~=0,j)>0) || ...
     all(model_ir.S(model_ir.S(:,j)~=0,j)<0)), 1:n));
importE    = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)>0), exCols));
exportE    = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)<0), exCols));
boundaryE  = unique([importE, exportE]);

% Initialize for rescue
covered      = false(n,1);
badPair_count = 0;

%% 

% Perform rescue
[allEFMs_missing, allX_missing, covered_rescue, badPair_rescue] = rescueUncoveredEFMs( ...
    template, aeq, Aineq, beq, bineq, ...
    lb0, ub0, vtype0, ...
    model_ir, covered, badPair_count, ...
    eps_flux, pPool, boundaryE );

% 6.4) Concatenate original and rescued results
allEFMs = [ allEFMs_original,   allEFMs_missing ];
allX    = [ allX,               allX_missing ];  % allX from enumeration is still in workspace

%% 

% 6.5) Build the full flux matrix
numEFMs_total = numel(allEFMs);
fluxMat = zeros(n, numEFMs_total);
for k = 1:numEFMs_total
    fluxMat(:,k) = allX{k};
end

% 6.6) Save combined results
save('EFMs_pool_badpair_afterPoolRescue_final.mat', ...
     'allEFMs', 'fluxMat', 'covered_rescue', 'badPair_rescue');

% Write combined CSV
Tfinal = array2table(fluxMat, ...
    'RowNames',      model_ir.rxns, ...
    'VariableNames', strcat("EFM", string(1:numEFMs_total)) );
writetable(Tfinal, 'EFMs_flux_final.csv', 'WriteRowNames', true);

fprintf('Rescue added %d EFMs. Total EFMs = %d. Results saved.\n', ...
        numel(allEFMs_missing), numEFMs_total);


%% === local helper functions ===

function tf = hasBadPair(support, rxnNames)
% true if any “_f”/“_b” pair appears together
    tf = false;
    for i=support
        rn = rxnNames{i};
        if endsWith(rn,'_f') && any(strcmp(rxnNames(support), [rn(1:end-2) '_b']))
            tf = true; return;
        end
    end
end

function template = applyNoGoodCut(template, supp, n)
% bans exactly the binary pattern in supp
    oth = setdiff(1:n, supp);
    cut = sparse(1,2*n);
    cut(n+supp) = 1;
    cut(n+oth)  = -1;
    template.A     = [template.A; cut];
    template.rhs   = [template.rhs; numel(supp)-1];
    template.sense = [template.sense; '<'];
end

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

function [allEFMs_missing, allX_missing, covered, badPair_count] = rescueUncoveredEFMs( ...
        template, aeq, Aineq, beq, bineq, lb0, ub0, vtype0, ...
        model_ir, covered, badPair_count, eps_flux, pPool, boundaryE )

    % boundaryE: vector of boundary reaction indices; pass [] to drop boundary req.
    m = size(aeq,1);
    n = numel(lb0)/2;

    % rebuild core inequality system (no implicit import/export cuts)
    Acore = Aineq;  bcore = bineq;
    if ~isempty(boundaryE)
        cutB         = sparse(1,2*n);
        cutB(n+boundaryE) = -1;
        Acore   = [Acore; cutB];
        bcore   = [bcore; -1];
    end

    % assemble static template
    template.A        = [aeq; Acore];
    template.rhs      = [beq; bcore];
    template.sense    = [ repmat('=',m,1); repmat('<',size(Acore,1),1) ];
    template.lb       = lb0;
    template.ub       = ub0;
    template.vtype    = vtype0;
    template.modelsense = 'min';
    template.obj      = [ zeros(n,1); ones(n,1) ];

    % solver params: single‐shot, silent
    params = pPool;
    params.PoolSearchMode = 0;
    params.PoolSolutions  = 1;
    params.OutputFlag     = 0;
    params.LogToConsole   = 0;
    params.LogFile        = '';

    stillUncovered  = [1015,1962];    % or pass in as additional argument
    allEFMs_missing = {};
    allX_missing    = {};

    for uc = stillUncovered(:).'
        fprintf('\n=== Rescue for reaction %d (%s) ===\n', uc, model_ir.rxns{uc});
        T = template;
        T.lb(n+uc) = 1;   % force reaction uc on

        noGoodA = []; noGoodb = [];
        found   = false;

        for attempt = 1:500
            % append any no‐good cuts
            if ~isempty(noGoodA)
                T.A     = [T.A; noGoodA];
                T.rhs   = [T.rhs; noGoodb];
                T.sense = [T.sense; repmat('<',size(noGoodA,1),1)];
            end

            sol = gurobi(T, params);
            if ~isfield(sol,'status') || ~strcmp(sol.status,'OPTIMAL')
                fprintf('  no more candidates (status=%s)\n', sol.status);
                break;
            end

            x    = sol.x;
            v    = x(1:n);
            z    = x(n+1:end);
            supp = find(z > eps_flux);

            % nullity via null(...) for consistent tol
            Ns      = null(full(model_ir.S(:,supp)));
            nullDim = size(Ns,2);

            badPairF = hasBadPair(supp, model_ir.rxns);

            fprintf('  cand |supp|=%3d, nullDim=%d, badPair=%d\n', ...
                    numel(supp), nullDim, badPairF);

            if nullDim==1 && ~badPairF
                fprintf('    → accepted EFM #%d (|supp|=%d)\n', ...
                        numel(allEFMs_missing)+1, numel(supp));
                allEFMs_missing{end+1} = supp;
                allX_missing   {end+1} = v;       % record the flux vector
                covered(supp)         = true;
                badPair_count         = badPair_count + badPairF;
                found = true;
                break;
            end

            % ban this exact support
            cut = sparse(1,2*n);
            cut(n+supp) = 1;
            cut(n+setdiff(1:n,supp)) = -1;
            noGoodA = [noGoodA; cut];
            noGoodb = [noGoodb; numel(supp)-1];
        end

        if ~found
            warning('  Could not rescue reaction %d\n', uc);
        end
    end
end

%%
