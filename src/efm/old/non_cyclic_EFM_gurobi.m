%% 0) Determine number of threads for intlinprog (OpenMP inside Gurobi)
nThreads = feature('numcores') - 1;
fprintf('Detected %d CPU threads; using multi-threaded Gurobi\n', nThreads);
setenv('OMP_NUM_THREADS', num2str(nThreads));

gurobiFolder = fullfile(getenv('GUROBI_HOME'),'examples','matlab');
addpath(gurobiFolder,'-begin');
fprintf('Calling Gurobi intlinprog wrapper from:\n  %s\n', which('intlinprog'));

%% 1) COBRA toolbox init
cobraModelFile      = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/iML1515.mat';
restrictToSubsystem = false;
subsystemName       = 'Glycolysis/Gluconeogenesis';
maxModes            = 5000;
eps_flux            = 1e-3;
M                   = 1e3;

%% 2) Load & optionally slice subsystem
orig = readCbModel(cobraModelFile);
if restrictToSubsystem
    model = sliceModel(orig, subsystemName);
else
    model = orig;
end

%% 3) Convert to irreversible & get sizes
model_ir = convertToIrreversible(model);
[m, n]    = size(model_ir.S);
nOrig     = numel(model.rxns);

%% 4) Build MILP matrices (static)
intcon = (n+1):(2*n);
aeq    = [model_ir.S, zeros(m,n)];  beq = zeros(m,1);
Aineq = [ eye(n),      M*eye(n);
         -eye(n),   -eps_flux*eye(n);
          zeros(1,n), ones(1,n) ];
bineq = [ M*ones(n,1);
         -eps_flux*ones(n,1);
          n-2 ];
lb = [model_ir.lb; zeros(n,1)];
ub = [model_ir.ub; ones(n,1)];

%% 5) Static constraints
% 5.1) split‐reaction exclusivity
revIdx = find(orig.lb<0 & orig.ub>0);
for k = 1:numel(revIdx)
    cut = zeros(1,2*n);
    cut(n + revIdx(k))     = -1;
    cut(n + nOrig + k)     = -1;
    Aineq = [Aineq; cut];
    bineq = [bineq; -1];
end
% 5.2 & 5.3) import/export & cycle‐prevention
exCols  = find(arrayfun(@(j) any(model_ir.S(:,j)~=0) && ...
               (all(model_ir.S(model_ir.S(:,j)~=0,j)>0) || ...
                all(model_ir.S(model_ir.S(:,j)~=0,j)<0)), 1:n));
importE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)>0), exCols));
exportE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)<0), exCols));
cutImp = zeros(1,2*n); cutImp(n+importE)=1;
cutExp = zeros(1,2*n); cutExp(n+exportE)=1;
Aineq  = [Aineq; cutImp; cutExp];
bineq  = [bineq; numel(importE)-7; numel(exportE)-7];

%% 6) Configure Gurobi‐via‐intlinprog options
opts = struct('Display','off');
opts.CutGeneration       = 'advanced';
opts.Heuristics          = 0.8;   
opts.MaxTime             = 3000;
opts.MaxNodes            = 1e6;
opts.RelativeGapTolerance= 5e-2;

%% 7) Prep progress & coverage tracking
hWait     = waitbar(0,'Enumerating EFMs...');
tStart    = tic;
efms      = cell(1, maxModes);
EFVs      = zeros(n, maxModes);
modeCount = 0;
covered   = false(n,1);          % ← track “seen” reactions

% *keep* the original objective: minimize total support
f = [ zeros(n,1);
      ones(n,1) ];
useCount = zeros(n,1);
%% 8) Enumerate EFMs, stop when all covered
while modeCount < maxModes
    [x,~,exitflag] = intlinprog(-f, intcon, Aineq, bineq, aeq, beq, lb, ub, [], opts);
    if exitflag ~= 1
        fprintf('Solver exitflag = %d; stopping.\n', exitflag);
        break;
    end
%after we found a support we restrict the search space to reactions which
%have not yet been seleted. This can be for example just int he space of
%export and impirt reactions which will invitably cause other pathways in
%the network to be triggered. 
    v    = x(1:n);
    supp = find(v >= eps_flux);

    % nullity check
    if numel(supp) - rank(model_ir.S(:, supp)) ~= 1
        warning('Invalid EFM—skipping.');
        cut = zeros(1,2*n); cut(n+supp)=-1;
        Aineq = [Aineq; cut]; bineq = [bineq; -1];
        continue;
    end

    % record & update coverage
    modeCount = modeCount + 1;
    efms{modeCount}      = supp;
    EFVs(supp,modeCount) = round(v(supp),6);
    covered(supp)        = true;
    useCount(supp) = useCount(supp) + 1;
    f = [ zeros(n,1);
          ones(n,1) + useCount ];
    % progress
    elapsed = toc(tStart);
    fprintf('[%5.1f s] Mode %d (|supp|=%d), covered %d/%d\n', ...
            elapsed, modeCount, numel(supp), sum(covered), n);

    if all(covered)
        fprintf('► All %d reactions covered after %d EFMs.\n', n, modeCount);
        break;
    end

    % exclusion cut
    cut = zeros(1,2*n); cut(n+supp)=-1;
    Aineq = [Aineq; cut]; bineq = [bineq; -1];

    waitbar(modeCount/maxModes, hWait, sprintf('Found %d of %d',modeCount,maxModes));
end

close(hWait);
fprintf('Done: %d EFMs in %.1f s.\n', modeCount, toc(tStart));

%% 9) Save to CSV
outDir  = '/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/';
outFile = fullfile(outDir,'efm_fluxes_coverage.csv');
fid = fopen(outFile,'w');
fprintf(fid,'Reaction');
for k=1:modeCount, fprintf(fid,',Mode%d',k); end
fprintf(fid,'\n');
for i=1:n
    fprintf(fid,'%s', model_ir.rxns{i});
    for k=1:modeCount
        fprintf(fid,',%.6f', EFVs(i,k));
    end
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('CSV saved: %s\n', outFile);
