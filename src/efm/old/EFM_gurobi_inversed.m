% 0) Gurobi threads & paths
nThreads = feature('numcores') ;
fprintf('Detected %d CPU threads; using multi-threaded Gurobi\n', nThreads);
setenv('OMP_NUM_THREADS', num2str(nThreads));
gurobiFolder = fullfile(getenv('GUROBI_HOME'),'examples','matlab');
addpath(gurobiFolder,'-begin');
fprintf('Calling Gurobi intlinprog wrapper from:\n  %s\n', which('intlinprog'));

% 1) Load the *irreversible* pruned model
data    = load("/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/pruned_iML1515_irrev_split_biomassrm_v4.mat",'pruned_ir');
model_ir = data.pruned_ir;   % already irreversible
[m, n] = size(model_ir.S);
%the second block of variables z (formerly y) will now mean z=1 active, z=0 inactive

% 3) Initialize enumeration parameters
eps_flux = 9e-6;
M        = 1e4;
EFVs     = zeros(n,n);

% 4) Build the static MILP matrices
aeq = [ model_ir.S, zeros(m,n) ];
beq = zeros(m,1);

%coupling: if z_i=0 -> v_i=0; if z_i=1 -> eps <= v_i <= M
Aineq = [ ...
    eye(n),      -M*eye(n);      %  v - M·z ≤ 0   → v ≤ M·z
   -eye(n),       eps_flux*eye(n)% -v + ε·z ≤ 0   → v ≥ ε·z
    %zeros(1,n),  -ones(1,n)       % -sum(z) ≤ -2   → sum(z) ≥ 2 (at least two actives)
];
bineq = [ ...
    zeros(n,1)
    zeros(n,1)
   %-2
];

lb0    = [ model_ir.lb; zeros(n,1) ];   % v ≥ LB,  z ≥ 0
ub0    = [ model_ir.ub; ones(n,1) ];    % v ≤ UB,  z ≤ 1
lb0(1067)=1;
intcon = (n+1):(2*n);

%% 5) Static constraints
%%5.1) split‐reaction exclusivity: z_f + z_b ≤ 1
%for i = 1:n
%    rxnName = model_ir.rxns{i};
%    if endsWith(rxnName, '_f')
%        baseName = rxnName(1:end-2);
%        bwdName  = [baseName, '_b'];
%        j = find(strcmp(model_ir.rxns, bwdName), 1);
%        if ~isempty(j)
%            cut = zeros(1, 2*n);
%            cut(n + i) = 1;    % +z_f
%            cut(n + j) = 1;    % +z_b
%            Aineq = [Aineq; cut];
%            bineq = [bineq; 1];
%        end
%    end
%end
%
%5.2 & 5.3) import/export & cycle‐prevention 
exCols  = find(arrayfun(@(j) any(model_ir.S(:,j)~=0) && ...
           (all(model_ir.S(model_ir.S(:,j)~=0,j)>0) || ...
            all(model_ir.S(model_ir.S(:,j)~=0,j)<0)), 1:n));
importE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)>0), exCols));

disp(numel(importE))
exportE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)<0), exCols));
%at least one import active: sum(z(importE)) ≥ 1 → -sum(z(importE)) ≤ -1
disp(numel(exportE))
cutImp  = zeros(1,2*n);
cutImp(n+importE) = -1;
%at least one export active: sum(z(exportE)) ≥ 1 → -sum(z(exportE)) ≤ -1
cutExp  = zeros(1,2*n);
cutExp(n+exportE) = -1;
Aineq   = [Aineq; cutImp; cutExp];
bineq   = [bineq; -1; -1];

% 6) Configure Gurobi‐via‐intlinprog options
opts = optimoptions('intlinprog', ...
    'Display',              'off',   ...  % no screen output
    'MaxTime',              3000,    ...  % 3000 s time limit
    'ConstraintTolerance',1e-9, ...
    'IntegerTolerance',   1e-9, ...
    'MaxFeasiblePoints',    Inf,     ...  % no limit on feasible solutions
    'RelativeGapTolerance', 5e-1,    ... 
    'AbsoluteGapTolerance', 0        ...  
);

% 7) Prep progress & coverage tracking
hWait     = waitbar(0,'Enumerating EFMs...');
tStart    = tic;
efms      = cell(n,1);
covered   = false(n,1);
f         = [ zeros(n,1); ones(n,1) ];  % minimize sum(z) → minimal active set
maxRetries_nullity=50;
maxRetries=3;
badpair_count=0;
% 8) Loop: force each reaction on (z_r=1), find one EFM
for r = 1:n
    if all(covered)
        fprintf('All reactions covered—terminating enumeration.\n');
        break;
    end
    if covered(r)
        fprintf('Skipping reaction %d (already covered).\n', r);
        continue;
    end

    fprintf('Computing EFM for reaction %d/%d ...\n', r, n);

    %Local copies so we can adjust without touching global
    eps_local = eps_flux;
    gap_local = opts.RelativeGapTolerance;

    %Save static constraints for resetting on each retry
    Aineq_static = Aineq;
    bineq_static = bineq;

    %Initialize bounds
    Acur = Aineq_static;
    bcur = bineq_static;
    lb = lb0;
    ub = ub0;
    %force reaction r active: z_r = 1
    lb(n + r) = 1;
    %ensure v_r ≥ eps
    if r~=1067
        lb(r)= eps_local;
    else 
        lb(r)=lb(r);
    end 
    
    retry  = 0;
    solved = false;
    nullity_retry=0;
    badPair = false;
    badPair_already=false;
    while ~solved && retry < maxRetries
        [x, ~, exitflag, output] = intlinprog(f, intcon, Acur, bcur, aeq, beq, lb, ub, [], opts);

        if exitflag == 1 && nullity_retry <= maxRetries_nullity
            v = x(1:n);
            z = x(n+1:2*n);
            badPair=false;
            %Extract support from z
            supp = find(z == 1);
            suppRxns = model_ir.rxns(supp);
%
            %% Quick scan for any forward/backward pair in the support
            
            badPair = false;

            % Scan for any forward/backward pair in the support
            for k = 1:numel(suppRxns)
                rxn = suppRxns{k};
                if endsWith(rxn, '_f')
                    base    = rxn(1:end-2);
                    partner = [base, '_b'];
                    % check if the backward partner is also in suppRxns
                    if any(strcmp(suppRxns, partner))
                        badPair = true;
                        break;
                    end
                end
            end
            
            if badPair
                badpair_count = badpair_count + 1;
                fprintf('    Detected forward/backward pair in support. Excluding only this exact support and retrying (badPair #%d).\n', badpair_count);
            
                others = setdiff(1:n, supp);
            
                % build the no-good cut:  sum(z(supp)) - sum(z(others)) <= |supp| - 1
                cut = zeros(1, 2*n);
                cut(n + supp)   =  1;   % +z_i for i in the support
                cut(n + others) = -1;   % -z_j for j outside the support
            
                Acur = [Acur; cut];
                bcur = [bcur; numel(supp) - 1];
            
            
                continue;
            end
            if numel(supp) - rank(model_ir.S(:, supp)) == 1
                %Valid EFM
                efms{r}       = supp;
                flux_r        = v(r);
                
                if any(abs(v(supp)) < 2e-6)
                    fprintf('One reaction has not enough flux\n');
                    %cut      = zeros(1,2*n);
                    %cut(n + supp) = 1;  % exclude these z's
                    %Acur    = [Acur; cut];
                    %bcur    = [bcur; numel(supp)-1];
                    others = setdiff(1:n, supp);
                    cut      = zeros(1,2*n);
                    cut(n + supp)  =  1;    % +z_i for i in supp
                    cut(n + others)= -1;    % –z_j for all j not in supp
                    Acur   = [Acur; cut];
                    bcur  = [bcur; numel(supp)-1];
                    nullity_retry=nullity_retry+1;
                    continue;
                end 
                %for k = 1:numel(suppRxns)
                %    rxn = suppRxns{k};
                %    if endsWith(rxn,'_f')
                %        base = rxn(1:end-2);
                %        if any(strcmp(suppRxns, [base '_b']))
                %            badPair = true;  
                %        end
                %    end
                %end
                %if badPair && ~badPair_already
                %    fprintf('We got a bad pair\n');
                %    badpair_count=badpair_count+1;
                %    badPair_already=true;
                %end
                fprintf('Flux through reaction %d (%s) = %.6g. This is valid.\n', r, model_ir.rxns{r}, flux_r);
                covered(supp) = true;
                solved        = true;
                EFVs(:,r)     = v;
                nullity_retry = 0;
                %cut      = zeros(1,2*n);
                %cut(n + supp) = 1;  % exclude these z's
                %Acur    = [Acur; cut];
                %bcur    = [bcur; numel(supp)-1];
                others = setdiff(1:n, supp);
                cut      = zeros(1,2*n);
                cut(n + supp)  =  1;    % +z_i for i in supp
                cut(n + others)= -1;    % –z_j for all j not in supp
                Acur   = [Acur; cut];
                bcur  = [bcur; numel(supp)-1];
                break;
            else
               
                %Nullity check failed → exclude this support and all
                %supersets
                others = setdiff(1:n, supp);
                cut      = zeros(1,2*n);
                cut(n + supp)  =  1;    % +z_i for i in supp
                cut(n + others)= -1;    % –z_j for all j not in supp
                Acur   = [Acur; cut];
                bcur  = [bcur; numel(supp)-1];
                continue;
            end
        end

        %exitflag ~= 1: infeasible or timed‐out → adjust and retry
        retry = retry + 1;
        nullity_retry=0;
        fprintf('  → intlinprog exitflag=%d (%s). Attempt %d/%d...\n', exitflag, output.message, retry, maxRetries);

        if retry == 1 
            %First retry: increase gap tolerance
            gap_local = 5e-1;
            opts.RelativeGapTolerance = gap_local;
            fprintf('    Increasing RelativeGapTolerance to %.2e and retrying...\n', gap_local);
            continue;

        elseif retry == 2
            %Second retry: lower eps_flux
            %eps_local = 1e-6;
            %Acur(1:n, n+1:2*n) = eps_local * eye(n);
            %lb(r)     = eps_local;
            %fprintf('    Lowering eps_flux to %.2e and retrying...\n', eps_local);
            %continue;
            gap_local = 9e-1;
            opts.RelativeGapTolerance = gap_local;
            fprintf('    Increasing RelativeGapTolerance to %.2e and retrying...\n', gap_local);
            continue;
        else
            warning('Failed to find EFM for reaction %d after %d retries. Excluding support.\n', r, retry);
            covered(r) = false;
            break;
        end
    end

    if solved
        %Restore gap and eps for next r
        opts.RelativeGapTolerance = 5e-2;
    end

    elapsed = toc(tStart);
    fprintf('[%5.1f s] Reaction %d done (|supp|=%d), covered %d/%d\n', elapsed, r, numel(efms{r}), sum(covered), n);
    waitbar(sum(covered)/n, hWait);
end

close(hWait);

%%

% 9) Collect reactions not covered
notCovered = find(~covered);
nUncov     = numel(notCovered);
fprintf('Uncovered reactions (%d):\n', numel(notCovered));
for k = 1:numel(notCovered)
    idx = notCovered(k);
    fprintf('%4d: %s\n', idx, model_ir.rxns{idx});
end
%%
coveredUncov = false(nUncov,1);
for idx = 1:nUncov
    if all(coveredUncov)
        fprintf('All reactions covered—terminating enumeration.\n');
        break;
    end
    if coveredUncov(idx)
        fprintf('Skipping reaction %d (already covered).\n', idx);
        continue;
    end
    r = notCovered(idx);
    fprintf('Computing EFM for reaction %d/%d ...\n', r, n);
    
    %Local copies so we can adjust without touching global
    eps_local = eps_flux;
    gap_local = opts.RelativeGapTolerance;

    %Save static constraints for resetting on each retry
    Aineq_static = Aineq;
    bineq_static = bineq;

    %Initialize bounds
    Acur = Aineq_static;
    bcur = bineq_static;
    lb = lb0;
    ub = ub0;
    %force reaction r active: z_r = 1
    lb(n + r) = 1;    %ensure v_r ≥ ep
    if r~=1067
        lb(r)= eps_local;
    else 
        lb(r)=lb(r);
    end 
    
    retry  = 0;
    solved = false;
    nullity_retry=0;
    badPair = false;
    while ~solved && retry < maxRetries
        [x, ~, exitflag, output] = intlinprog(f, intcon, Acur, bcur, aeq, beq, lb, ub, [], opts);

        if exitflag == 1 && nullity_retry <= maxRetries_nullity
            v = x(1:n);
            z = x(n+1:2*n);
            badPair=false;
            %Extract support from z
            supp = find(z == 1);
            suppRxns = model_ir.rxns(supp);
%
            %% Quick scan for any forward/backward pair in the support
            badPair = false;

            % Scan for any forward/backward pair in the support
            for k = 1:numel(suppRxns)
                rxn = suppRxns{k};
                if endsWith(rxn, '_f')
                    base    = rxn(1:end-2);
                    partner = [base, '_b'];
                    % check if the backward partner is also in suppRxns
                    if any(strcmp(suppRxns, partner))
                        badPair = true;
                        break;
                    end
                end
            end
            
            if badPair
                badpair_count = badpair_count + 1;
                fprintf('    Detected forward/backward pair in support. Excluding only this exact support and retrying (badPair #%d).\n', badpair_count);
            
                others = setdiff(1:n, supp);
            
                % build the no-good cut:  sum(z(supp)) - sum(z(others)) <= |supp| - 1
                cut = zeros(1, 2*n);
                cut(n + supp)   =  1;   % +z_i for i in the support
                cut(n + others) = -1;   % -z_j for j outside the support
            
                Acur = [Acur; cut];
                bcur = [bcur; numel(supp) - 1];
            
            
                continue;
            end
            if numel(supp) - rank(model_ir.S(:, supp)) == 1
                %Valid EFM
                efms{r}       = supp;
                flux_r        = v(r);
                
                if any(abs(v(supp)) < 2e-6)
                    fprintf('One reaction has not enough flux\n');
                    %cut      = zeros(1,2*n);
                    %cut(n + supp) = 1;  % exclude these z's
                    %Acur    = [Acur; cut];
                    %bcur    = [bcur; numel(supp)-1];
                    others = setdiff(1:n, supp);
                    cut      = zeros(1,2*n);
                    cut(n + supp)  =  1;    % +z_i for i in supp
                    cut(n + others)= -1;    % –z_j for all j not in supp
                    Acur   = [Acur; cut];
                    bcur  = [bcur; numel(supp)-1];
                    nullity_retry=nullity_retry+1;
                    continue;
                end 
                %for k = 1:numel(suppRxns)
                %    rxn = suppRxns{k};
                %    if endsWith(rxn,'_f')
                %        base = rxn(1:end-2);
                %        if any(strcmp(suppRxns, [base '_b']))
                %            badPair = true;  
                %        end
                %    end
                %end
                %if badPair
                %    fprintf('We got a bad pair\n');
                %    badpair_count=badpair_count+1;
                %end
                fprintf('Flux through reaction %d (%s) = %.6g. This is valid.\n', r, model_ir.rxns{r}, flux_r);
                covered(supp) = true;
                coveredUncov(idx) = true;
                solved        = true;
                EFVs(:,r)     = v;
                nullity_retry = 0;
                %cut      = zeros(1,2*n);
                %cut(n + supp) = 1;  % exclude these z's
                %Acur    = [Acur; cut];
                %bcur    = [bcur; numel(supp)-1];
                others = setdiff(1:n, supp);
                cut      = zeros(1,2*n);
                cut(n + supp)  =  1;    % +z_i for i in supp
                cut(n + others)= -1;    % –z_j for all j not in supp
                Acur   = [Acur; cut];
                bcur  = [bcur; numel(supp)-1];
                break;
            else
                %Nullity check failed → exclude support and retry
                %cut      = zeros(1,2*n);
                %cut(n + supp) = 1;  % exclude these z's
                %Acur    = [Acur; cut];
                %bcur    = [bcur; numel(supp)-1];
                others = setdiff(1:n, supp);
                cut      = zeros(1,2*n);
                cut(n + supp)  =  1;    % +z_i for i in supp
                cut(n + others)= -1;    % –z_j for all j not in supp
                Acur   = [Acur; cut];
                bcur  = [bcur; numel(supp)-1];
                nullity_retry=nullity_retry+1;
                continue;
            end
        end
        %exitflag ~= 1: infeasible or timed‐out → adjust and retry
        retry = retry + 1;
        nullity_retry=0;
        fprintf('  → intlinprog exitflag=%d (%s). Attempt %d/%d...\n', exitflag, output.message, retry, maxRetries);

        if retry == 1 
            %First retry: increase gap tolerance
            gap_local = 5e-1;
            opts.RelativeGapTolerance = gap_local;
            fprintf('    Increasing RelativeGapTolerance to %.2e and retrying...\n', gap_local);
            continue;

        elseif retry == 2
            %Second retry: lower eps_flux
            %eps_local = 1e-6;
            %Acur(1:n, n+1:2*n) = eps_local * eye(n);
            %lb(r)     = eps_local;
            %fprintf('    Lowering eps_flux to %.2e and retrying...\n', eps_local);
            %continue;
            gap_local = 9e-1;
            opts.RelativeGapTolerance = gap_local;
            fprintf('    Increasing RelativeGapTolerance to %.2e and retrying...\n', gap_local);
            continue;
        else
            warning('Failed to find EFM for reaction %d after %d retries. Excluding support.\n', r, retry);
            covered(r) = false;
            break;
        end
    end

    if solved
        %Restore gap and eps for next r
        opts.RelativeGapTolerance = 5e-2;
    end

    
    fprintf(' Reaction %d done (|supp|=%d), covered %d/%d\n', r, numel(efms{r}), sum(covered), n);
    

end 

%%
stillUncov = notCovered(~coveredUncov);
fprintf('Still uncovered after two passes (%d):\n', numel(stillUncov));
for r = stillUncov(:)'
    fprintf('  %4d: %s\n', r, model_ir.rxns{r});
end

%%

save( ...
  '/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/EFMs_by_reaction_v3_inverse.mat', ...
  'efms', 'covered', 'stillUncov', 'EFVs' );
fprintf('All results saved to EFMs_by_reaction.mat\n');
fprintf('Enumeration complete: %d/%d reactions covered.\n', sum(~cellfun(@isempty,efms)), n);

%% 12) Save CSV of reaction‐flux matrix
outDir  = '/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/ecoli/';
outFile = fullfile(outDir, 'efm_reaction_fluxes_2_inverse.csv');
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