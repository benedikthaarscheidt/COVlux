%% 0) Determine number of threads for intlinprog
nThreads = feature('numcores');
fprintf('Detected %d CPU threads; using multi-threaded built-in intlinprog.\n', nThreads);
setenv('OMP_NUM_THREADS', num2str(nThreads));
maxNumCompThreads(1);  % prevent MATLAB BLAS from oversubscribing

%% 1) COBRA toolbox init
cobraModelFile      = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/iML1515.mat';
restrictToSubsystem = false;
subsystemName       = 'Glycolysis/Gluconeogenesis';
maxModes            = 3;
eps_flux            = 1e-3; % small threshold
M                   = 1e3;



%% 2) Load & optionally slice subsystem
orig = readCbModel(cobraModelFile);

% Exchange-reaction counts
nRegexEx   = sum(~cellfun(@isempty, regexp(orig.rxns,'^EX_')));
Sorig      = orig.S;
n0         = size(Sorig,2);
stoichMask = arrayfun(@(j) any(Sorig(:,j)~=0) && ...
               (all(Sorig(Sorig(:,j)~=0,j)>0) || all(Sorig(Sorig(:,j)~=0,j)<0)), ...
               1:n0);
nStoichEx  = sum(stoichMask);
nRevStoich = sum(orig.lb(stoichMask)<0 & orig.ub(stoichMask)>0);

fprintf('Regex-based:   %d exchange rxns\n', nRegexEx);
fprintf('Stoich-based:  %d exchange rxns (%d reversible)\n', nStoichEx, nRevStoich);

if restrictToSubsystem
    model = sliceModel(orig, subsystemName);
    fprintf('Sliced to "%s" (%d rxns)\n', subsystemName, numel(model.rxns));
else
    model = orig;
    fprintf('Using full model (%d rxns)\n', numel(model.rxns));
end

%% 3) Convert to irreversible & get sizes
model_ir = convertToIrreversible(model);
[m, n]    = size(model_ir.S);
nOrig     = numel(model.rxns);
fprintf('Irreversible: %d reactions (%d orig + %d split), %d metabolites\n', ...
        n, nOrig, n-nOrig, m);

%% 4) Build MILP matrices
intcon = (n+1):(2*n);
f      = [zeros(n,1); ones(n,1)];  % minimize support

aeq = [model_ir.S, zeros(m,n)];  beq = zeros(m,1);
Aineq = [ eye(n),      M*eye(n);
         -eye(n),   -eps_flux*eye(n);
          zeros(1,n), ones(1,n) ];
bineq = [ M*ones(n,1);
         -eps_flux*ones(n,1);
          n-2 ];
lb = [model_ir.lb; zeros(n,1)];
ub = [model_ir.ub; ones(n,1)];

%% 5) Static constraints
revIdx = find(orig.lb<0 & orig.ub>0);
for k=1:numel(revIdx)
    cut = zeros(1,2*n);
    cut(n + revIdx(k))     = -1;
    cut(n + nOrig + k)     = -1;
    Aineq = [Aineq; cut];
    bineq = [bineq; -1];
end

exCols  = find(arrayfun(@(j) any(model_ir.S(:,j)~=0) && ...
               (all(model_ir.S(model_ir.S(:,j)~=0,j)>0) || ...
                all(model_ir.S(model_ir.S(:,j)~=0,j)<0)), 1:n));
importE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)>0), exCols));
exportE = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)<0), exCols));
cutImp  = zeros(1,2*n); cutImp(n+importE)=1;
cutExp  = zeros(1,2*n); cutExp(n+exportE)=1;
Aineq   = [Aineq; cutImp; cutExp];
bineq   = [bineq; numel(importE)-1; numel(exportE)-1];

% Locate your Gurobi wrapper folder (as set in startup.m)
gurobiFolder = fullfile(getenv('GUROBI_HOME'),'examples','matlab');

% If the Gurobi wrapper is ahead on the path, remove it temporarily
wrapperPath = which('intlinprog');
if contains(wrapperPath, gurobiFolder)
    fprintf('Removing Gurobi wrapper from path to call built-in intlinprog.\n');
    rmpath(gurobiFolder);
end

% Verify we’re about to use the built-in version
wrapperPath = which('intlinprog');
fprintf('Calling built-in intlinprog from:\n  %s\n', wrapperPath);

% Only use valid built-in options
optsBI = optimoptions('intlinprog', ...
    'Display',              'off', ...
    'MaxNodes',             1e5, ...
    'RelativeGapTolerance', 1e-6);




%% 7) Prepare progress bar & storage
hWait     = waitbar(0,'Enumerating EFMs...');
tStart    = tic;
efms      = cell(1,maxModes);
EFVs      = zeros(n,maxModes);
modeCount = 0;

%% 8) Enumerate EFMs with dynamic exclusion using built-in solver
while modeCount < maxModes
    % Force call to built-in intlinprog, bypassing any wrapper
    [x,~,exitflag] = builtin('intlinprog', ...
       -f, intcon, Aineq, bineq, aeq, beq, lb, ub, opts);

    if exitflag ~= 1
        fprintf('Solver exitflag = %d; stopping enumeration.\n', exitflag);
        break;
    end

    v    = x(1:n);
    supp = find(v >= eps_flux);

    % verify nullity
    nullity = numel(supp) - rank(model_ir.S(:,supp));
    if nullity ~= 1
        warning('Invalid EFM (|supp|=%d): nullity=%d — skipping.', numel(supp), nullity);
        cut = zeros(1,2*n); cut(n+supp)=-1;
        Aineq = [Aineq; cut]; bineq=[bineq;-1];
        continue;
    end

    modeCount = modeCount + 1;
    efms{modeCount}      = supp;
    EFVs(supp,modeCount) = round(v(supp),6);

    fprintf('[%5.1f s] Mode %d/%d (|supp|=%d)\n', ...
            toc(tStart), modeCount, maxModes, numel(supp));

    cut = zeros(1,2*n); cut(n+supp)=-1;
    Aineq = [Aineq; cut]; bineq=[bineq;-1];

    waitbar(modeCount/maxModes, hWait);
end

close(hWait);
fprintf('Done: %d EFMs in %.1f s.\n', modeCount, toc(tStart));

%% 9) Save to CSV
outputFile = fullfile('/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/', ...
                      'efm_fluxes_output_builtin.csv');
fid = fopen(outputFile,'w');
fprintf(fid,'Reaction');
for k = 1:modeCount, fprintf(fid,',Mode%d', k); end
fprintf(fid,'\n');
for i = 1:n
    fprintf(fid,'%s',model_ir.rxns{i});
    for k = 1:modeCount
        fprintf(fid,',%.6f',EFVs(i,k));
    end
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('CSV saved: %s\n', outputFile);
