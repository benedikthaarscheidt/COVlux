% ------------------------------------------------------------------------
% MILP‐based minimal‐support modes on Recon3D
%   • optional subsystem restriction
%   • enforce v=0 when y=1, v≥ε when y=0
%   • require ≥2 active reactions
%   • progress bar via waitbar()
% ------------------------------------------------------------------------

%% ——— USER SETUP ———————————————————————————————
cobraModelFile      = 'Recon3D_301.mat';
restrictToSubsystem = true;                         % false = full model
subsystemName       = 'Glycolysis/Gluconeogenesis';
maxModes            = 20;                           % how many modes to extract at max
eps_flux            = 1e-3;                          % min flux when active
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
    if ~any(mask), error('No reactions in "%s".', subsystemName); end
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
[m,n]   = size(model_ir.S);
intcon  = (n+1):(2*n);              % y‐variables
f       = [zeros(n,1); ones(n,1)];  % maximize ∑y → minimize #active

% (a) Steady‐state: S·v = 0
Aeq = [model_ir.S, zeros(m,n)];
beq = zeros(m,1);

% (b) Link v & y: v + M·y ≤ M
Aineq = [eye(n),   M*eye(n)];
bineq = M * ones(n,1);

% (c) Enforce v ≥ ε·(1−y):  −v − ε·y ≤ −ε
Aineq = [Aineq;
         -eye(n), -eps_flux*eye(n)];
bineq = [bineq;
         -eps_flux*ones(n,1)];

% (d) At least two active reactions: ∑y ≤ n−2
Aineq = [Aineq;
         zeros(1,n), ones(1,n)];
bineq = [bineq;
         n-2];

% Bounds
lb = [model_ir.lb; zeros(n,1)];
ub = [model_ir.ub; ones(n,1)];

%% 4) Enumerate minimal‐support candidates with progress bar
efms = {};
opts = optimoptions('intlinprog','Display','off');

hWait = waitbar(0, 'Enumerating EFMs...');

while numel(efms) < maxModes
    [x,~,exitflag] = intlinprog(-f,intcon,Aineq,bineq,Aeq,beq,lb,ub,opts);
    if exitflag ~= 1
        break;
    end

    y    = round(x(n+1:end));
    supp = find(y==0);
    efms{end+1} = supp;

    % update progress
    waitbar(numel(efms)/maxModes, hWait, ...
        sprintf('Found %d of %d EFMs...', numel(efms), maxModes));

    % forbid this exact support: −∑_{i∈supp} y_i ≤ −1
    cut = zeros(1,2*n);
    cut(n+supp) = -1;
    Aineq = [Aineq; cut];
    bineq = [bineq; -1];
end

close(hWait);

%% 5) Print the resulting modes
if restrictToSubsystem
    fprintf('\nFound %d minimal‐support modes in "%s":\n', numel(efms), subsystemName);
else
    fprintf('\nFound %d minimal‐support modes in full model:\n', numel(efms));
end
for k = 1:numel(efms)
    rxnList = model_ir.rxns(efms{k});
    fprintf('  Mode %2d (|supp|=%d): %s\n', k, numel(efms{k}), strjoin(rxnList, ', '));
end
