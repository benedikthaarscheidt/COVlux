%% 0) Gurobi threads & paths
nThreads = feature('numcores') - 1;
fprintf('Detected %d CPU threads; using multi-threaded Gurobi\n', nThreads);
setenv('OMP_NUM_THREADS', num2str(nThreads));
gurobiFolder = fullfile(getenv('GUROBI_HOME'),'examples','matlab');
addpath(gurobiFolder,'-begin');
fprintf('Calling Gurobi intlinprog wrapper from:\n  %s\n', which('intlinprog'));

%% 1) Load the *irreversible* pruned model
data    = load('/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/pruned_iML1515_irrev_split_biomassrm.mat','pruned_ir');
model_ir = data.pruned_ir;   % already irreversibl

n = numel(model_ir.rxns);
[m, n] = size(model_ir.S);

%% 3) Initialize enumeration parameters
eps_flux = 1e-5;
M        = 1e6;
EFVs     = zeros(n,n);

%% 4) Build the static MILP matrices
% 4.1) S·v = 0
aeq = [ model_ir.S, zeros(m,n) ];
beq = zeros(m,1);

% 4.2) big-M flux‐activation + cardinality ≥ 2 active
Aineq = [ ...
    eye(n),        M*eye(n);      %  v + M·y ≤ M
   -eye(n),  -eps_flux*eye(n);    % −v − ε·y ≤ −ε
    zeros(1,n),    ones(1,n)      % sum(y) ≤ n−2
];
bineq = [ ...
    M*ones(n,1);
   -eps_flux*ones(n,1);
    n-2
];

% 4.3) bounds and integer indices
lb0    = [ model_ir.lb; zeros(n,1) ];
ub0    = [ model_ir.ub; ones(n,1) ];
intcon = (n+1):(2*n);

%% 5) Static constraints
% 5.1) split‐reaction exclusivity (match “_f” ↔ “_b”)
for i = 1:n
    rxnName = model_ir.rxns{i};
    if endsWith(rxnName, '_f')
        baseName = rxnName(1:end-2);
        bwdName  = [baseName, '_b'];
        j = find(strcmp(model_ir.rxns, bwdName), 1);
        if ~isempty(j)
            cut = zeros(1, 2*n);
            cut(n + i) = -1;   % –y_fwd
            cut(n + j) = -1;   % –y_bwd
            Aineq = [Aineq; cut];
            bineq = [bineq; -1];
        end
    end
end

% 5.2 & 5.3) import/export & cycle‐prevention (unchanged)
exCols  = find(arrayfun(@(j) any(model_ir.S(:,j)~=0) && ...
           (all(model_ir.S(model_ir.S(:,j)~=0,j)>0) || ...
            all(model_ir.S(model_ir.S(:,j)~=0,j)<0)), 1:n));
importE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)>0), exCols));
exportE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)<0), exCols));
cutImp  = zeros(1,2*n);
cutImp(n+importE)=1;
cutExp  = zeros(1,2*n);
cutExp(n+exportE)=1;
Aineq   = [Aineq; cutImp; cutExp];
bineq   = [bineq; numel(importE)-1; numel(exportE)-1];

%% 6) Configure Gurobi‐via‐intlinprog options
opts = struct('Display','off');
opts.CutGeneration        = 'advanced';
opts.Heuristics           = 0.8;   
opts.MaxTime              = 3000;
opts.MaxNodes             = 1e6;
opts.RelativeGapTolerance = 5e-2;
opts.FeasibilityTol = 1e-8;

%% 7) Prep progress & coverage tracking
hWait     = waitbar(0,'Enumerating EFMs...');
tStart    = tic;
efms      = cell(n,1);
covered   = false(n,1);
f         = [ zeros(n,1); ones(n,1) ];  % MAX sum(y) = #inactives

%% 8) Loop: force each reaction on (y_r=0), find one EFM
for r = 1:n
    if covered(r)
        fprintf('Skipping reaction %d (already covered).\n', r);
        continue;
    end

    fprintf('Computing EFM for reaction %d/%d …\n', r, n);

    % reset MILP
    Acur = Aineq;
    bcur = bineq;
    lb   = lb0;
    ub   = ub0;

    % force reaction r active (y_r = 0)
    ub(n + r) = 0;
    lb(r)     = eps_flux;
    % try up to 3 times to solve the MILP
    maxRetries = 3;
    retry      = 0;
    solved     = false;

    while true
        [x, ~, exitflag, output] = intlinprog(-f, intcon, Acur, bcur, aeq, beq, lb, ub, [], opts);

        if exitflag == 1
            % we’ve got a solution—check support
            v    = x(1:n);
            y = x(n+1:2*n);

            % Extract support from y
            supp = find(y == 0);

            if numel(supp) - rank(model_ir.S(:, supp)) == 1
                % valid EFM
                efms{r}       = supp;
                covered(supp) = true;
                flux_r        = v(r);
                fprintf('Flux through reaction %d (%s) = %.6g\n',r, model_ir.rxns{r}, flux_r);
                covered(r)    = true;
                solved        = true;
                EFVs(:,r) = v; 
                break;
            else
                % not a true EFM—exclude this support and retry
                cut      = zeros(1,2*n);
                cut(n + supp) = -1;
                Acur    = [Acur; cut];
                bcur    = [bcur; -1];
                % continue without counting as a retry
                continue;
            end
        end

        % if we reach here, exitflag ~= 1 → a failure/time‐out
        retry = retry + 1;
        fprintf('  → intlinprog exitflag=%d (%s). Retry %d/%d…\n', ...
                exitflag, output.message, retry, maxRetries);

        if retry >= maxRetries
            warning('Still no solution for reaction %d after %d retries. Computing later with relaxed constraints', r, retry);
            covered(r) = false;
            break;
        end
        % otherwise loop and try again
    end

    % if we found a valid EFM, globally exclude it so it won’t repeat
    if solved
        cut              = zeros(1,2*n);
        cut(n + efms{r}) = -1;
        Aineq            = [Aineq; cut];
        bineq            = [bineq; -1];
    end

    % progress
    elapsed = toc(tStart);
    fprintf('[%5.1f s] Reaction %d done (|supp|=%d), covered %d/%d\n', ...
            elapsed, r, numel(efms{r}), sum(covered), n);
    waitbar(sum(covered)/n, hWait);
end

close(hWait);

%% 9) Collect reactions not covered
notCovered = find(~covered);
fprintf('Reactions not covered after initial pass (%d):\n', numel(notCovered));
disp(model_ir.rxns(notCovered));

%% 10) Second pass: try uncovered reactions with relaxed constraints
% Identify rows of the original cardinality & import/export cuts:
rowCard = 2*n + 1;       % the "sum(y) ≤ n-2" row
rowImp  = rowCard + 1;   % import cut
rowExp  = rowCard + 2;   % export cut

% Build relaxed constraint set (drop card, imp & exp)
A_relaxed = Aineq;
b_relaxed = bineq;
rowsToDrop = sort([rowCard, rowImp, rowExp],'descend');
A_relaxed(rowsToDrop,:) = [];
b_relaxed(rowsToDrop)   = [];

hWait2 = waitbar(0,'Relaxed enumeration for uncovered reactions...');
for k = 1:numel(notCovered)
    r = notCovered(k);
    if covered(r)
        continue;
    end

    fprintf('Relaxed try for reaction %d …\n', r);
    % reset MILP
    Acur = A_relaxed;
    bcur = b_relaxed;
    lb   = lb0;
    ub   = ub0;

    ub(n + r) = 0;  % force on

    % single attempt (you can also wrap retries if desired)
    [x, ~, exitflag] = intlinprog(-f, intcon, Acur, bcur, aeq, beq, lb, ub, [], opts);
    if exitflag == 1
        v    = x(1:n);
        supp = find(v >= eps_flux);
        if numel(supp) - rank(model_ir.S(:, supp)) == 1
            efms{r}       = supp;
            covered(supp) = true;
            fprintf('  → Found relaxed EFM for reaction %d (|supp|=%d).\n', r, numel(supp));
            % exclude globally to avoid repeats
            cut              = zeros(1,2*n);
            cut(n + supp)    = -1;
            Aineq            = [Aineq; cut];
            bineq            = [bineq; -1];
        else
            fprintf('  → Relaxed solution for reaction %d failed nullity check.\n', r);
        end
    else
        fprintf('  → No EFM found under relaxed constraints for reaction %d (exitflag=%d).\n', r, exitflag);
    end

    waitbar(k/numel(notCovered), hWait2);
end
close(hWait2);

%% 11) Final report and save
finalNotCovered = find(~covered);
fprintf('\nFinal uncovered reactions (%d):\n', numel(finalNotCovered));
disp(model_ir.rxns(finalNotCovered));

save('/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/EFMs_by_reaction.mat', ...
     'efms', 'covered', 'notCovered', 'finalNotCovered');
fprintf('All results saved to EFMs_by_reaction.mat\n');
fprintf('Enumeration complete: %d/%d reactions covered.\n', sum(~cellfun(@isempty,efms)), n);

%% 9) Save results
outDir  = '/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/ecoli/';
outFile = fullfile(outDir,'efm_reaction_fluxes_1.csv');
fid     = fopen(outFile,'w');

% Header
fprintf(fid,'Reaction');
for k = 1:n
    fprintf(fid,',Mode%d', k);
end
fprintf(fid,'\n');

% Rows
for i = 1:n
    fprintf(fid,'%s', model_ir.rxns{i});
    for k = 1:n
        fprintf(fid,',%.6g', EFVs(i,k));
    end
    fprintf(fid,'\n');
end

fclose(fid);
fprintf('CSV saved: %s\n', outFile);
