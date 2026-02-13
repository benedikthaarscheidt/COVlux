%% 0) Gurobi threads & paths
nThreads = feature('numcores') - 1;
fprintf('Detected %d CPU threads; using multi-threaded Gurobi\n', nThreads);
setenv('OMP_NUM_THREADS', num2str(nThreads));
gurobiFolder = fullfile(getenv('GUROBI_HOME'),'examples','matlab');
addpath(gurobiFolder,'-begin');
fprintf('Calling Gurobi intlinprog wrapper from:\n  %s\n', which('intlinprog'));

%% 1) Load the *irreversible* pruned model
data    = load("/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/pruned_iML1515_irrev_split_biomassrm.mat",'pruned_ir');
model_ir = data.pruned_ir;   % already irreversible
[m, n] = size(model_ir.S);
%model_ir.lb(1067)=0;
%% 3) Initialize enumeration parameters  ← (renumbered now that n has changed)
eps_flux = 9e-6;
M        = 9e6;
EFVs     = zeros(n,n);

%% 4) Build the static MILP matrices
% (exactly as before, but now n includes ‘DM_succoa’)
aeq = [ model_ir.S, zeros(m,n) ];
beq = zeros(m,1);

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

lb0    = [ model_ir.lb; zeros(n,1) ];
ub0    = [ model_ir.ub; ones(n,1) ];
ub0(n+1067)=0;
intcon = (n+1):(2*n);
%% 5) Static constraints
% 5.1) split‐reaction exclusivity (match “_f” ↔ “_b”)
%for i = 1:n
%    rxnName = model_ir.rxns{i};
%    if endsWith(rxnName, '_f')
%        baseName = rxnName(1:end-2);
%        bwdName  = [baseName, '_b'];
%        j = find(strcmp(model_ir.rxns, bwdName), 1);
%        if ~isempty(j)
%            cut = zeros(1, 2*n);
%            cut(n + i) = -1;   % –y_fwd
%            cut(n + j) = -1;   % –y_bwd
%            Aineq = [Aineq; cut];
%            bineq = [bineq; -1];
%        end
%    end
%end

% 5.2 & 5.3) import/export & cycle‐prevention 
exCols  = find(arrayfun(@(j) any(model_ir.S(:,j)~=0) && ...
           (all(model_ir.S(model_ir.S(:,j)~=0,j)>0) || ...
            all(model_ir.S(model_ir.S(:,j)~=0,j)<0)), 1:n));
importE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)>0), exCols));
disp(importE)
exportE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)<0), exCols));
cutImp  = zeros(1,2*n);
cutImp(n+importE)=1;
cutExp  = zeros(1,2*n);
cutExp(n+exportE)=1;
Aineq   = [Aineq; cutImp; cutExp];
bineq   = [bineq; numel(importE)-1; numel(exportE)-1];

%% 6) Configure Gurobi‐via‐intlinprog options
opts = optimoptions('intlinprog', ...
    'Display',              'off',   ...  % no screen output
    'MaxTime',              3000,    ...  % 3000 s time limit
    'MaxFeasiblePoints',    inf,     ...  % no limit on feasible solutions
    'RelativeGapTolerance', 5e-2,    ...  % 5% MIP‐gap
    'AbsoluteGapTolerance', 0        ...  
);
%% 7) Prep progress & coverage tracking
hWait     = waitbar(0,'Enumerating EFMs...');
tStart    = tic;
efms      = cell(n,1);
covered   = false(n,1);
f         = [ zeros(n,1); ones(n,1) ];  % MAX sum(y) = inactives
maxRetries_nullity=50;
maxRetries=3;
badpair_count=0;
%% 8) Loop: force each reaction on (y_r=0), find one EFM
for r = 1:n
    if all(covered)
        fprintf('All reactions covered—terminating enumeration.\n');
        break;
    end
    if covered(r)
        fprintf('Skipping reaction %d (already covered).\n', r);
        continue;
    end

    fprintf('Computing EFM for reaction %d/%d …\n', r, n);

    % Local copies so we can adjust without touching global
    eps_local = eps_flux;                     % first-attempt threshold
    gap_local = opts.RelativeGapTolerance;

    % Save static constraints for resetting on each retry
    Aineq_static = Aineq;
    bineq_static = bineq;

    % First‐try: eps_local = 9e-6
    Acur = Aineq_static;
    bcur = bineq_static;
    %Acur(end, r)   = 1;          % v_r ≥ eps_local  →  −v_r ≤ −eps_local
    %Acur(end, n+r) = eps_local;
    %bcur(end)      = eps_local;
    lb = lb0; ub = ub0;
    ub(n + r) = 0;
    if r~=1067
        lb(r)= eps_local;
    else 
        lb(r)=lb0(r);
    end 

    retry  = 0;
    solved = false;
    lastSupp = [];  % remember last tested support-
    nullity_retry=0;
    
    while ~solved && retry < maxRetries
        [x, ~, exitflag, output] = intlinprog(-f, intcon, Acur, bcur, aeq, beq, lb, ub, [], opts);

        if exitflag == 1 && nullity_retry <= maxRetries_nullity;
            v = x(1:n);
            y = x(n+1:2*n);
            
            % Extract support from y
            supp = find(y == 0);
            suppRxns = model_ir.rxns(supp);
            
            for k = 1:numel(suppRxns)
                rxn = suppRxns{k};
                if endsWith(rxn,'_f')
                    base = rxn(1:end-2);
                    if any(strcmp(suppRxns, [base '_b']))
                        badpair_count=badpair_count+1;
                        badPair = true;

                        fprintf('We got a bad pair\n');
                    end
                end
            end



            %% Quick scan for any forward/backward pair in the support
            %badPair = false;
            %for k = 1:numel(suppRxns)
            %    rxn = suppRxns{k};
            %    if endsWith(rxn,'_f')
            %        base = rxn(1:end-2);
            %        if any(strcmp(suppRxns, [base '_b']))
            %            badPair = true;
            %        end
            %    end
            %end
            %
            %if badPair
            %    % Exclude this support and retry
            %    cut       = zeros(1,2*n);
            %    cut(n + supp) = -1;       % sum_{i∈supp} y_i ≤ |supp|−1
            %    Acur      = [Acur; cut];
            %    bcur      = [bcur; -1];
            %    nullity_retry = nullity_retry + 1;
            %    continue
            %end

            if numel(supp) - rank(model_ir.S(:, supp)) == 1
                % Valid EFM
                efms{r}       = supp;
                
                flux_r        = v(r);
                fprintf('Flux through reaction %d (%s) = %.6g\n', r, model_ir.rxns{r}, flux_r);
                if flux_r<eps_flux
                    fprintf('fugging flux too small. we try again. \n');
                    cut      = zeros(1,2*n);
                    cut(n + supp) = -1;
                    Acur    = [Acur; cut];
                    bcur    = [bcur; -1];
                    nullity_retry=nullity_retry+1;
                    continue;
                end 

                
                %fprintf('Support reactions (indices and names):\n');
                %for i = 1:numel(supp)
                %    idx = supp(i);
                %    fprintf('  %3d: %s\n', idx, model_ir.rxns{idx});
                %end
                covered(supp) = true;
                
                solved        = true;
                EFVs(:,r)     = v;
                nullity_retry = 0;
                break;
            else
                %fprintf('    Nullity check failed. Retrying (%d/%d)…\n', nullity_retry, maxRetries_nullity);
                cut      = zeros(1,2*n);
                cut(n + supp) = -1;
                Acur    = [Acur; cut];
                bcur    = [bcur; -1];
                nullity_retry=nullity_retry+1;
                continue;
            end
        end

        % exitflag ~= 1: infeasible or timed‐out → adjust and retry
        retry = retry + 1;
        nullity_retry=0;
        fprintf('  → intlinprog exitflag=%d (%s). Attempt %d/%d…\n', ...
                exitflag, output.message, retry, maxRetries);

        if retry == 1 
            % First retry: increase gap tolerance
            gap_local = 5e-1;
            opts.RelativeGapTolerance = gap_local;
            fprintf('    Increasing RelativeGapTolerance to %.2e and retrying…\n', gap_local);
            %Acur = [Aineq_static; zeros(1, 2*n)];
            %bcur = [bineq_static; 0];
            %Acur(end, r)   = 1;
            %Acur(end, n+r) = eps_local;
            %bcur(end)      = eps_local;
            %lb(r)          = eps_local;
            %ub(n+r)        = 0;
            continue;

        elseif retry == 2
            % Second retry: lower eps_flux
            %eps_local = 1e-6;
            %Acur(n+1:2*n, n+1:2*n) = -eps_local * eye(n);
            %bcur(n+1:2*n)=-eps_local;
            %fprintf('    Lowering eps_flux to %.2e and retrying…\n', eps_local);
            %Acur = [Aineq_static; zeros(1, 2*n)];
            %bcur = [bineq_static; 0];
            %Acur(end, r)   = 1;
            %Acur(end, n+r) = eps_local;
            %bcur(end)      = eps_local;
            %lb(r)          = eps_local;
            %ub(n+r)        = 0;
            gap_local = 9e-1;
            opts.RelativeGapTolerance = gap_local;
            fprintf('    Increasing RelativeGapTolerance to %.2e and retrying…\n', gap_local);
            continue;

        else
            % Third retry: give up on r, now exclude last failed support
            warning('Failed to find EFM for reaction %d after %d retries. Excluding support.\n', r, retry);
            covered(r) = false;
            break;
        end
    end

    if solved
        % Exclude this EFM globally
        cut             = zeros(1, 2*n);
        cut(n + efms{r}) = -1;
        Aineq           = [Aineq; cut];
        bineq           = [bineq; -1];
        % Restore gap and eps for next r
        opts.RelativeGapTolerance = gap_local;
        eps_local = eps_flux;
    end

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
rowFluxMinus_rel = (n+1) : (2*n);

% Build the original A_relaxed_static / b_relaxed_static at eps_flux=9e-6
R       = 0;  % no split‐reaction cuts in irreversible model
rowCard = 2*n + 1;       
rowImp  = rowCard + 1;   
rowExp  = rowCard + 2;   

A_relaxed_static = Aineq;
b_relaxed_static = bineq;
rowsToDrop       = sort([rowCard, rowImp, rowExp], 'descend');
A_relaxed_static(rowsToDrop, :) = [];
b_relaxed_static(rowsToDrop)   = [];

hWait2 = waitbar(0,'Relaxed enumeration for uncovered reactions...');
for k = 1:numel(notCovered)
    r = notCovered(k);
    if covered(r)
        continue;
    end

    fprintf('Relaxed try for reaction %d …\n', r);

    % Initialize local eps and gap for this reaction
    eps_local = eps_flux;                   % 9e-6 on first attempt
    gap_local = opts.RelativeGapTolerance;  % 5e-4 as set above

    % Build a fresh copy of A_relaxed_static and b_relaxed_static 
    % which currently use eps_flux = 9e-6 in their “–v – eps·y ≤ –eps” rows.
    Arel_static = A_relaxed_static;
    brel_static = b_relaxed_static;

    % Ensure the “–v – eps·y ≤ –eps” block is correct for eps_local = 9e-6:
    Arel_static(rowFluxMinus_rel, 1:n      ) = -eye(n);
    Arel_static(rowFluxMinus_rel, n+1:2*n  ) = -eps_local * eye(n);
    brel_static(rowFluxMinus_rel)            = -eps_local * ones(n,1);

    retry = 0;
    nullity_retry = 0;
    solved = false;

    while ~solved && retry < maxRetries
        % (Re)build Acur / bcur from the static relaxed block
        Acur = [Arel_static; zeros(1, 2*n)];
        bcur = [brel_static;    0      ];

        % Append “force reaction r on” row: v_r + eps_local·y_r ≥ eps_local
        Acur(end,   r    ) =  1;
        Acur(end, n + r  ) =  eps_local;
        bcur(end      ) =  eps_local;

        % Reset bounds
        lb = lb0;
        ub = ub0;
        ub(n + r) = 0;      % y_r = 0  → reaction r active
        lb(r)     = eps_local;  % v_r ≥ eps_local

        % Solve
        [x, ~, exitflag, output] = intlinprog(-f, intcon, Acur, bcur, aeq, beq, lb, ub, [], opts);

        if exitflag == 1 && nullity_retry <= maxRetries_nullity
            v    = x(1:n);
            supp = find(v >= eps_local);

            if numel(supp) - rank(model_ir.S(:, supp)) == 1
                % Found a valid EFM under relaxed constraints
                efms{r}       = supp;
                covered(supp) = true;
                fprintf('  → Found relaxed EFM for reaction %d (|supp|=%d).\n', r, numel(supp));
                % Exclude it globally so future loops skip it
                cut           = zeros(1, 2*n);
                cut(n + supp) = -1;
                Aineq         = [Aineq; cut];
                bineq         = [bineq; -1];
                solved = true;
                break;

            else
                % Nullity check failed → exclude support and retry (no retry++ for gap)
                fprintf('    Relaxed nullity check failed. Retrying (exclusion) …\n');
                cut = zeros(1, 2*n);
                cut(n + supp) = -1;
                Arel_static = [Arel_static; cut];
                brel_static = [brel_static; -1];
                nullity_retry = nullity_retry + 1;
                continue;
            end
        end

        % exitflag ≠ 1: infeasible or timed‐out
        retry = retry + 1;
        nullity_retry = 0;
        fprintf('  → Relaxed intlinprog exitflag=%d (%s). Attempt %d/%d…\n', ...
                exitflag, output.message, retry, maxRetries);

        if retry == 1
            % First retry: loosen gap tolerance, keep eps_local = 9e-6
            gap_local = 5e-2;
            opts.RelativeGapTolerance = gap_local;
            fprintf('    Increasing RelativeGapTolerance to %.2e and retrying…\n', gap_local);
            continue;

        elseif retry == 2
            % Second retry: lower eps_local to 1e-6, update static block
            eps_local = 1e-6;
            lb(r)     = eps_local;
            ub(n+r)   = 0;

            % Overwrite static “–v – eps·y ≤ –eps” rows:
            Arel_static(rowFluxMinus_rel, 1:n      ) = -eye(n);
            Arel_static(rowFluxMinus_rel, n+1:2*n  ) = -eps_local * eye(n);
            brel_static(rowFluxMinus_rel)            = -eps_local * ones(n,1);

            fprintf('    Lowering eps_flux to %.2e and retrying…\n', eps_local);
            continue;

        else
            % Third retry: give up on this reaction
            warning('  → Giving up on relaxed enumeration for reaction %d after %d retries.\n', r, retry);
            covered(r) = false;
            break;
        end
    end
    waitbar(k/numel(notCovered), hWait2);
end
close(hWait2);

%% 11) Final report and save
finalNotCovered = find(~covered);
fprintf('\nFinal uncovered reactions (%d):\n', numel(finalNotCovered));
disp(model_ir.rxns(finalNotCovered));

save('/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/EFMs_by_reaction_v2.mat', ...
     'efms', 'covered', 'notCovered', 'finalNotCovered');
fprintf('All results saved to EFMs_by_reaction.mat\n');
fprintf('Enumeration complete: %d/%d reactions covered.\n', sum(~cellfun(@isempty,efms)), n);

%% 12) Save CSV of reaction‐flux matrix
outDir  = '/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/ecoli/';
outFile = fullfile(outDir, 'efm_reaction_fluxes_1.csv');
fid     = fopen(outFile, 'w');

% Header
fprintf(fid, 'Reaction');
for k = 1:n
    fprintf(fid, ',Mode%d', k);
end
fprintf(fid, '\n');

% Rows
for i = 1:n
    fprintf(fid, '%s', model_ir.rxns{i});
    for k = 1:n
        fprintf(fid, ',%.6g', EFVs(i, k));
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('CSV saved: %s\n', outFile);
