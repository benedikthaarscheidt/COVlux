%% 0) Initialize COBRA toolbox
fprintf('Starting COBRA toolbox init...\n');
cobraPath = '/work/haarscheid/master_thesis/cobratoolbox';
addpath(genpath(cobraPath));
initCobraToolbox(false);

%% 1) Determine number of workers from SLURM (and clamp to cluster max)
requestedCores = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(requestedCores) || requestedCores < 1
    requestedCores = feature('numcores');
end

% use the Slurm profile for both clamping and pool-launching
pctCluster = parcluster('slurm');      
maxWorkers  = pctCluster.NumWorkers;    % cluster‐wide max
nCores      = min(requestedCores, maxWorkers);

if requestedCores > maxWorkers
    fprintf('Warning: requested %d cores, but Slurm profile max is %d; using %d cores\n', ...
            requestedCores, maxWorkers, nCores);
else
    fprintf('Launching with %d workers\n', nCores);
end

%% 2) Start parallel pool via Slurm
pool = gcp('nocreate');
if isempty(pool) || pool.NumWorkers ~= nCores
    if ~isempty(pool)
        delete(pool);
    end
    pool = parpool(pctCluster, nCores);   % submits workers into your SLURM allocation
end



%% 3) Load & optionally slice subsystem

cobraModelFile      = '/work/haarscheid/master_thesis/Models/generic_models/E_coli/iML1515.mat';
restrictToSubsystem = false;
subsystemName       = 'Glycolysis/Gluconeogenesis';
maxModes            = 30;
eps_flux            = 1e-3; % small threshold
M                   = 1e3;

orig = readCbModel(cobraModelFile);

% --- EXCHANGE counts in original model
% 1) regex count
regexIdx = ~cellfun(@isempty, regexp(orig.rxns, '^EX_'));
nRegexEx = sum(regexIdx);
fprintf('Regex-based: found %d exchange reactions in original model.', nRegexEx);

% 2) stoichiometry-based count --> only verificaction
Sorig = orig.S;
n0    = size(Sorig,2);
stoichMask = arrayfun(@(j) any(Sorig(:,j)~=0) && ...
    (all(Sorig(Sorig(:,j)~=0,j)>0) || all(Sorig(Sorig(:,j)~=0,j)<0)), 1:n0);
stoichIdx = find(stoichMask);
fprintf('Stoich-based: found %d exchange reactions in original model.', numel(stoichIdx));

% 3) among these, check reversibility via lb and ub --> this is just verification 
nRevStoich = sum(orig.lb(stoichIdx)<0 & orig.ub(stoichIdx)>0);
fprintf('Of these, %d are reversible (lb<0 & ub>0).', nRevStoich);

if restrictToSubsystem
    model = sliceModel(orig, subsystemName);
    fprintf('Sliced to subsystem "%s" with %d reactions.\n', subsystemName, numel(model.rxns));
else
    model = orig;
    fprintf('Using full model with %d reactions.\n', numel(orig.rxns));
end

%% 3.2) Convert to irreversible and get sizes
model_ir = convertToIrreversible(model);
[m, n]    = size(model_ir.S);
nOrig     = numel(model.rxns);
fprintf('Converted to irreversible: %d reactions (%d original + %d split), %d metabolites.\n', ...
    n, nOrig, n-nOrig, m);

%% 4) Build MILP matrices
intcon = (n+1):(2*n);
f      = [zeros(n,1); ones(n,1)];    % minimize support via sum(y)

aeq = [model_ir.S, zeros(m,n)];
beq = zeros(m,1);

Aineq = [ eye(n),      M*eye(n);
         -eye(n),   -eps_flux*eye(n);
          zeros(1,n), ones(1,n) ];
bineq = [ M*ones(n,1);
         -eps_flux*ones(n,1);
          n - 2           ];

lb = [model_ir.lb; zeros(n,1)];
ub = [model_ir.ub; ones(n,1)];   

%% 5) Static constraints
% 5.1) split-reaction exclusivity: forbid both directions ON for all
% reversible reactions as this would contradic thermodnamic feasibility 
revIdx = find(orig.lb < 0 & orig.ub > 0);
for k = 1:numel(revIdx)
    fwd = revIdx(k);
    bwd = nOrig + k;
    cut = zeros(1,2*n);
    cut(n + fwd) = -1;
    cut(n + bwd) = -1;
    Aineq = [Aineq; cut];
    bineq = [bineq; -1];
end

% 5.2) classify import/export exchange reactions
fprintf('Classifying import/export exchange reactions...\n');
exCols = find(arrayfun(@(j) any(model_ir.S(:,j)~=0) && ...
    (all(model_ir.S(model_ir.S(:,j)~=0,j)>0) || all(model_ir.S(model_ir.S(:,j)~=0,j)<0)), 1:n));
importEx = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)>0), exCols));
exportEx = exCols(arrayfun(@(j) all(model_ir.S(model_ir.S(:,j)~=0,j)<0), exCols));
fprintf('Found %d import and %d export exchange reactions.\n', numel(importEx), numel(exportEx));

% 5.3) at least one import AND one export per mode --> this prevents
% internal cycles 
cutImp = zeros(1,2*n);
cutImp(n + importEx) = 1;
Aineq = [Aineq; cutImp];
bineq = [bineq; numel(importEx)-1];

cutExp = zeros(1,2*n);
cutExp(n + exportEx) = 1;
Aineq = [Aineq; cutExp];
bineq = [bineq; numel(exportEx)-1];

%% 6) Configure intlinprog
opts = optimoptions('intlinprog', 'Display','off', ...
    'Heuristics','advanced','CutGeneration','advanced','MaxNodes',1e5);

%% 7) Prepare progress bar & storage
tStart = tic;
efms = {};
EFVs = zeros(n, maxModes);
modeCount = 0;

%% 8) Enumerate EFMs with dynamic exclusion
while modeCount < maxModes
    [x,~,exitflag] = intlinprog(-f, intcon, Aineq, bineq, aeq, beq, lb, ub, opts);
    if exitflag ~= 1
        fprintf('Solver exitflag = %d; stopping enumeration.\n', exitflag);
        break;
    end
    v    = x(1:n);
    supp = find(v >= eps_flux);

    % verify nullity==1 to ensure we actually get valid modes
    subS    = model_ir.S(:, supp);
    rnk     = rank(subS);
    nullity = numel(supp) - rnk;
    if nullity ~= 1
        warning('Invalid EFM candidate (support size %d): nullity = %d. Skipping.', numel(supp), nullity);
        cut = zeros(1,2*n);
        cut(n + supp) = -1;
        Aineq = [Aineq; cut];
        bineq = [bineq; -1];
        continue;
    end

    % store mode
    modeCount = modeCount + 1;
    efms{modeCount}       = supp;
    EFVs(supp, modeCount) = round(v(supp),6);
    names = model_ir.rxns(supp);
    fprintf('  Mode %d found (|supp| = %d) with reactions:\n', modeCount, numel(supp));
    fprintf('    %s\n', strjoin(names, ', '));

    % dynamic exclusion cut so that we do not pick the exact same EFM
    % again.
    cut = zeros(1,2*n);
    cut(n + supp) = -1;
    Aineq = [Aineq; cut];
    bineq = [bineq; -1];


end

totalTime = toc(tStart);
fprintf('Done: %d EFMs in %.1f s.\n', modeCount, totalTime);

%% 10) Save to CSV
outFile = fullfile(pwd, 'efm_fluxes_output.csv');
fid     = fopen(outFile, 'w');
fprintf(fid, 'Reaction');
for k = 1:found
    fprintf(fid, ',Mode%d', k);
end
fprintf(fid, '\n');
for i = 1:n
    fprintf(fid, '%s', model_ir.rxns{i});
    for k = 1:found
        fprintf(fid, ',%.6f', EFVs(i,k));
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('CSV saved: %s\n', outFile);

%% 11) Cleanup & exit
delete(pool);
exit(0);


