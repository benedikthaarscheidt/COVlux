% ------------------------------------------------------------------------
% MILP‐based minimal‐support modes on Recon3D
%   • optional subsystem restriction
%   • enforce v=0 when y=1, v≥ε when y=0
%   • require ≥2 active reactions
%   • cancelable progress bar (Stop button) that immediately halts solver
%   • save partial results on Stop
% ------------------------------------------------------------------------

%% ——— USER SETUP ———————————————————————————————
cobraModelFile      = 'Recon3D_301.mat';
restrictToSubsystem = false;                         % false = full model
subsystemName       = 'Glycolysis/Gluconeogenesis'; % only used if restrictToSubsystem
maxModes            = 5;                            % maximum modes to extract
eps_flux            = 1e-3;                         % minimal flux when active
M                   = 1e3;                          % big‐M for linking
%% ——————————————————————————————————————————————

%% 1) Load & optionally slice to subsystem
orig = readCbModel(cobraModelFile);
if restrictToSubsystem
    ss = orig.subSystems;
    if ischar(ss), ss = cellstr(ss);
    elseif iscell(ss) && all(cellfun(@iscell,ss))
        ss = cellfun(@(c) strjoin(c,'; '), ss,'Uni',false);
    end
    mask = contains(ss, subsystemName,'IgnoreCase',true);
    if ~any(mask)
        error('No reactions found for "%s".', subsystemName);
    end
    fprintf('Restricting to %d reactions in "%s".\n', nnz(mask), subsystemName);
    model.S    = orig.S(:,mask);
    model.rxns = orig.rxns(mask);
    model.lb   = orig.lb(mask);
    model.ub   = orig.ub(mask);
    if isfield(orig,'c'), model.c = orig.c(mask); end
else
    fprintf('Using full model with %d reactions.\n', numel(orig.rxns));
    model = orig;
end

%% 2) Convert to irreversible
model_ir = convertToIrreversible(model);

%% 3) Build MILP matrices
[m, n]    = size(model_ir.S);
intcon    = (n+1):(2*n);             % y‐variable indices
f         = [zeros(n,1); ones(n,1)]; % maximize ∑ y → minimize #inactive

% (a) Steady‐state: S·v = 0
Aeq = [model_ir.S, zeros(m,n)];
beq = zeros(m,1);

% (b) Linking: v + M·y ≤ M
Aineq = [eye(n), M*eye(n)];
bineq = M*ones(n,1);

% (c) Minimal‐flux: v ≥ ε·(1−y) → −v − ε·y ≤ −ε
Aineq = [Aineq;
         -eye(n), -eps_flux*eye(n)];
bineq = [bineq;
         -eps_flux*ones(n,1)];

% (d) At least two active reactions: ∑y ≤ n−2
Aineq = [Aineq;
         zeros(1,n), ones(1,n)];
bineq = [bineq;
         n-2];

lb = [model_ir.lb; zeros(n,1)];
ub = [model_ir.ub; ones(n,1)];

%% 4) Enumerate minimal‐support candidates with cancelable waitbar
efms = {};
EFVs = zeros(n, maxModes);

% Create a waitbar with a Stop button
hWait = waitbar(0, 'Enumerating EFMs...', ...
    'Name', 'EFM Enumeration', ...
    'CreateCancelBtn', 'setappdata(gcbf, ''stop'', 1)');
setappdata(hWait, 'stop', 0);

% Set up intlinprog options with an OutputFcn to catch the Stop click
opts = optimoptions('intlinprog', ...
    'Display', 'off', ...
    'OutputFcn', @stopButtonOutputFcn);

while numel(efms) < maxModes
    % Check if Stop was clicked before starting the solver
    if getappdata(hWait, 'stop')
        break;
    end

    [x,~,exitflag] = intlinprog(-f, intcon, Aineq, bineq, Aeq, beq, lb, ub, opts);

    % If the solver was stopped or failed, break immediately
    if exitflag ~= 1
        break;
    end

    % Extract flux vector and binary y-vector
    v = x(1:n);
    y = round(x(n+1:end));

    % Define support: reactions where y==0 (active) AND v≥eps_flux
    supp = find((y==0) & (abs(v) >= eps_flux));

    % Store support and fluxes
    efms{end+1}          = supp;
    EFVs(:, numel(efms)) = v;

    % Update progress
    waitbar(numel(efms)/maxModes, hWait, ...
        sprintf('Found %d of %d EFMs...', numel(efms), maxModes));

    % Forbid this exact y-support in subsequent iterations (use y==0 indices)
    cut = zeros(1, 2*n);
    cut(n + find(y==0)) = -1;  % -sum(y==0) ≤ -1
    Aineq = [Aineq; cut];
    bineq = [bineq; -1];
end

% Close waitbar if it still exists
if isgraphics(hWait), delete(hWait); end

% Determine how many modes were actually found
nFound        = numel(efms);
reactionNames = model_ir.rxns;
varNames      = arrayfun(@(k) sprintf('Mode%d',k), 1:nFound, 'Uni', false);

%% 5) Trim flux matrix to true supports and save results e.g., zero out any non-support entries in each mode
for k = 1:nFound
    mask_k = false(n,1);
    mask_k(efms{k}) = true;
    tmp = EFVs(:,k);
    tmp(~mask_k) = 0;
    EFVs(:,k) = tmp;
end

% Round to 6 decimal places
EFVs = round(EFVs, 6);

% Write a CSV
outputFile = '/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_1.csv';
fid = fopen(outputFile, 'w');

% Header line
fprintf(fid, 'Reaction');
for k = 1:nFound
    fprintf(fid, ',%s', varNames{k});
end
fprintf(fid, '\n');

% Data lines
for i = 1:numel(reactionNames)
    fprintf(fid, '%s', reactionNames{i});
    for k = 1:nFound
        fprintf(fid, ',%.6f', EFVs(i,k));
    end
    fprintf(fid, '\n');
end
fclose(fid);

%% 6) Print summary of active reactions (y==0) with flux values
if restrictToSubsystem
    fprintf('\nFound %d minimal‐support modes in "%s":\n', nFound, subsystemName);
else
    fprintf('\nFound %d minimal‐support modes in full model:\n', nFound);
end
for k = 1:nFound
    supp = efms{k};
    names = model_ir.rxns(supp);
    vals  = EFVs(supp, k);
    entries = cellfun(@(nm,fl) sprintf('%s(%.6f)', nm, fl), names, num2cell(vals), 'Uni', false);
    fprintf('  Mode %2d: %s\n', k, strjoin(entries, ', '));
end

%% ——— Local stop button OutputFcn —————————————————
function stop = stopButtonOutputFcn(~, optimValues, ~)
    stop = false;
    h = findall(0, 'Type', 'figure', 'Name', 'EFM Enumeration');
    if ~isempty(h) && getappdata(h, 'stop')
        stop = true;
        delete(h);
    end
end
