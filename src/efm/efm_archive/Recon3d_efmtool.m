% ------------------------------------------------------------------------
% EFM enumeration via reaction‐formulas on Recon3D (optional subsystem)
% ------------------------------------------------------------------------

%% ——— USER SETUP ———————————————————————————————
cobraModelFile        = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/Recon3D_301/Recon3D_301.mat';
restrictToSubsystem   = true;                          % true = slice to subsystem; false = use full model
subsystemName         = 'Glycolysis/Gluconeogenesis';  % only used if restrictToSubsystem==true
maxTimeSeconds        = 7200;                          % timeout for EFMtool
%% ——————————————————————————————————————————————

%% 1) Load the full Recon3D model
origModel = readCbModel(cobraModelFile);

if restrictToSubsystem
    % Normalize subSystems into cellstr
    ssField = origModel.subSystems;
    if ischar(ssField)
        ssList = cellstr(ssField);
    elseif iscell(ssField) && all(cellfun(@iscell, ssField))
        ssList = cellfun(@(c) strjoin(c,'; '), ssField,'Uni',false);
    else
        ssList = ssField;
    end

    % Show available subsystems (optional)
    uniqueSubs = unique(ssList);
    disp('Available subsystems (first 50):');
    disp(uniqueSubs(1:min(end,50)));

    % Pick out reactions in the desired subsystem
    mask = contains(ssList, subsystemName, 'IgnoreCase', true);
    if ~any(mask)
        error('No reactions found for subsystem "%s".', subsystemName);
    end
    rxnsSub = origModel.rxns(mask);
    fprintf('Restricting to %d reactions in subsystem "%s".\n', nnz(mask), subsystemName);

    % Build mini‐model
    mini.S    = origModel.S(:, mask);
    mini.mets = origModel.mets;
    mini.rxns = rxnsSub;
    mini.lb   = origModel.lb(mask);
    mini.ub   = origModel.ub(mask);
    if isfield(origModel,'c')
        mini.c = origModel.c(mask);
    end
else
    % Use full model as mini‐model
    fprintf('Running on the full model with %d reactions.\n', numel(origModel.rxns));
    mini = origModel;
end

%% 2) Convert to irreversible format
miniIrrev = convertToIrreversible(mini);

%% 3) Build reaction‐formula strings
S     = miniIrrev.S;
mets  = miniIrrev.mets;
rnames= miniIrrev.rxns;
nRxn  = size(S,2);

rformulas = cell(1,nRxn);
for i = 1:nRxn
    negIdx    = find(S(:,i)<0);
    posIdx    = find(S(:,i)>0);
    negCoeffs = -full(S(negIdx,i));
    posCoeffs =  full(S(posIdx,i));
    if isempty(negIdx)
        lhs = {' '};
    else
        lhs = arrayfun(@(k) sprintf('%g %s', negCoeffs(k), mets{negIdx(k)}), ...
                       (1:numel(negIdx))', 'Uni',false);
    end
    if isempty(posIdx)
        rhs = {' '};
    else
        rhs = arrayfun(@(k) sprintf('%g %s', posCoeffs(k), mets{posIdx(k)}), ...
                       (1:numel(posIdx))', 'Uni',false);
    end
    rformulas{i} = sprintf('%s --> %s', strjoin(lhs,' + '), strjoin(rhs,' + '));
end

%% 4) Ensure ./tmp exists for I/O
if ~exist('tmp','dir')
    mkdir('tmp');
end

%% 5) Configure EFMtool options
opts = CreateFluxModeOpts( ...
    'arithmetic',  'fractional', ...
    'compression', 'off',        ... % valid: 'off' or 'gzip'
    'tmpdir',      'tmp',        ...
    'parallelize', true,         ...
    'timeout',     maxTimeSeconds ...
);

%% 6) Run EFMtool via reaction‐formulas interface
fprintf('Starting formula‐based EFM enumeration on %d reactions…\n', nRxn);
try
    t0   = tic;
    mnet = CalculateFluxModes(rformulas, rnames, opts);
    elapsed = toc(t0);

    if isstruct(mnet) && isfield(mnet,'numEFMs')
        fprintf('→ Done: found %d EFMs in %.2f seconds\n', mnet.numEFMs, elapsed);
        save('efm_results.mat','mnet','opts','rnames');
    else
        warning('EFMtool returned unexpected output:');
        disp(mnet);
    end
catch ME
    fprintf('EFMtool failed: %s\n', ME.message);
    disp(getReport(ME));
    return;
end

%% 7) Collect and display reactions for each EFM
nEFMs        = size(mnet.efms, 2);
efmReactions = cell(nEFMs,1);

for k = 1:nEFMs
    nz = find(mnet.efms(:,k) ~= 0);
    efmReactions{k} = rnames(nz);
end

for k = 1:nEFMs
    fprintf('EFM %2d: %s\n', k, strjoin(efmReactions{k}, ', '));
end
