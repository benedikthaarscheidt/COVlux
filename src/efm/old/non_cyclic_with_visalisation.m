%% 0) Determine number of threads for intlinprog
nThreads = feature('numcores');
fprintf('Detected %d CPU threads; using multi-threaded intlinprog.\n', nThreads);
setenv('OMP_NUM_THREADS', num2str(nThreads));

%% 1) SETUP
cobraModelFile      = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/e_coli_core.mat';
restrictToSubsystem = false;
subsystemName       = 'Glycolysis/Gluconeogenesis';
maxModes            = 80;
eps_flux            = 1e-3; % small threshold
M                   = 1e3;

%% 2) Load &  slice subsystem
orig = readCbModel(cobraModelFile);

% --- EXCHANGE counts in original model
% 1) regex-based count
regexIdx = ~cellfun(@isempty, regexp(orig.rxns, '^EX_'));
nRegexEx = sum(regexIdx);
fprintf('Regex-based: found %d exchange reactions in original model.', nRegexEx);

% 2) stoichiometry-based count
Sorig = orig.S;
n0    = size(Sorig,2);
stoichMask = arrayfun(@(j) any(Sorig(:,j)~=0) && ...
    (all(Sorig(Sorig(:,j)~=0,j)>0) || all(Sorig(Sorig(:,j)~=0,j)<0)), 1:n0);
stoichIdx = find(stoichMask);
fprintf('Stoich-based: found %d exchange reactions in original model.', numel(stoichIdx));

% 3) among these, check reversibility via bounds
nRevStoich = sum(orig.lb(stoichIdx)<0 & orig.ub(stoichIdx)>0);
fprintf('Of these, %d are reversible (lb<0 & ub>0).', nRevStoich);

if restrictToSubsystem
    model = sliceModel(orig, subsystemName);
    fprintf('Sliced to subsystem "%s" with %d reactions.\n', subsystemName, numel(model.rxns));
else
    model = orig;
    fprintf('Using full model with %d reactions.\n', numel(orig.rxns));
end

%% 3) Convert to irreversible and get sizes
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
ub = [model_ir.ub; ones(n,1)];   % enforce y ∈ {0,1}

%% 5) Static constraints
% 5.1) split-reaction exclusivity: forbid both directions ON
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

% 5.3) at least one import AND one export per mode
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
hWait = waitbar(0, 'Enumerating EFMs...');
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

    % verify nullity
    subS    = model_ir.S(:, supp);
    rnk     = rank(subS);
    nullity = numel(supp) - rnk;
    disp(nullity)
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

    % dynamic exclusion cut
    cut = zeros(1,2*n);
    cut(n + supp) = -1;
    Aineq = [Aineq; cut];
    bineq = [bineq; -1];

    % update waitbar
    waitbar(modeCount/maxModes, hWait, sprintf('Found %d of %d EFMs', modeCount, maxModes));
end

close(hWait);
totalTime = toc(tStart);
fprintf('Done: %d EFMs in %.1f s.\n', modeCount, totalTime);

%% 9) Save to CSV
outputFile = fullfile('/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/', 'efm_fluxes_output_ecoli_core.csv');
fid = fopen(outputFile,'w');
fprintf(fid, 'Reaction');
for k = 1:modeCount, fprintf(fid, ',Mode%d', k); end
fprintf(fid, '\n');
for i = 1:n
    fprintf(fid, '%s', model_ir.rxns{i});
    for k = 1:modeCount
        fprintf(fid, ',%.6f', EFVs(i,k));
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('CSV saved: %s\n', outputFile);

%% 10) Done
fprintf('All done!');

%% 11) Visualize as metabolite–reaction graph (metabolites=nodes, reactions=edges)
% We connect substrates to products for each irreversible reaction
Smat = model_ir.S;
metCount = size(Smat,1);
reacCount = size(Smat,2);
from = [];
to   = [];
edgeReac = zeros(0,1);
for j = 1:reacCount
    subs = find(Smat(:,j) < 0);
    prods = find(Smat(:,j) > 0);
    for si = subs'
        for pi = prods'
            from(end+1) = si; %#ok<AGROW>
            to(end+1)   = pi; %#ok<AGROW>
            edgeReac(end+1) = j; %#ok<AGROW>
        end
    end
end

% Build directed graph: nodes are metabolites, edges are reactions
MG = digraph(from, to, ones(size(from)), model_ir.mets);

% Plot full metabolic graph
figure;
h = plot(MG, 'Layout','layered', 'NodeLabel', MG.Nodes.Name);
title('Metabolite–Reaction Graph (all reactions)');
h.LineWidth = 0.5;

% Highlight edges for each mode
for modeIdx = 1:modeCount
    edgesToHighlight = ismember(edgeReac, efms{modeIdx});
    % reset all edges to light gray
    highlight(h, 'Edges', 1:numedges(MG), 'EdgeColor', [0.8 0.8 0.8], 'LineWidth', 0.5);
    % highlight this mode's reaction edges in red
    idx = find(edgesToHighlight);
    highlight(h, 'Edges', idx, 'EdgeColor', [1 0 0], 'LineWidth', 2);
    title(sprintf('Mode %d support highlighted (edges)', modeIdx));
    pause(3);
end
