% ------------------------------------------------------------------------
% EFM enumeration via reaction‐formulas on Recon3D (optional subsystem)
%   • stops after timeout
%   • saves mnet, plus a CSV of reaction coefficients per EFM
% ------------------------------------------------------------------------

%% ——— USER SETUP ———————————————————————————————
cobraModelFile      = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/e_coli_core_splitandpruned.mat';
restrictToSubsystem = false;                         
subsystemName       = 'Glycolysis/Gluconeogenesis';  
maxTimeSeconds      = 5;                             
%% ——————————————————————————————————————————————

%% 1) Load & optionally slice
origModel = readCbModel(cobraModelFile);
if restrictToSubsystem
    ssField = origModel.subSystems;
    if ischar(ssField)
        ssList = cellstr(ssField);
    elseif iscell(ssField) && all(cellfun(@iscell, ssField))
        ssList = cellfun(@(c) strjoin(c,'; '), ssField,'Uni',false);
    else
        ssList = ssField;
    end
    mask = contains(ssList, subsystemName,'IgnoreCase',true);
    if ~any(mask)
        error('No reactions found for "%s".', subsystemName);
    end
    fprintf('Restricting to %d reactions in "%s".\n', nnz(mask), subsystemName);
    mini.S    = origModel.S(:,mask);
    mini.mets = origModel.mets;
    mini.rxns = origModel.rxns(mask);
    mini.lb   = origModel.lb(mask);
    mini.ub   = origModel.ub(mask);
    if isfield(origModel,'c'), mini.c = origModel.c(mask); end
else
    fprintf('Running full model with %d reactions.\n', numel(origModel.rxns));
    mini = origModel;
end

%% 2) Convert to irreversible
miniIrrev = convertToIrreversible(mini);
rnames    = miniIrrev.rxns;
nRxn      = numel(rnames);

%% 3) Build reaction‐formula strings
S    = miniIrrev.S;
mets = miniIrrev.mets;
rformulas = cell(1,nRxn);
for i = 1:nRxn
    negIdx = find(S(:,i)<0);
    posIdx = find(S(:,i)>0);
    % left side
    if isempty(negIdx)
        lhsStr = ' ';
    else
        parts = cell(numel(negIdx),1);
        for j = 1:numel(negIdx)
            parts{j} = sprintf('%g %s', -full(S(negIdx(j),i)), mets{negIdx(j)});
        end
        lhsStr = strjoin(parts, ' + ');
    end
    % right side
    if isempty(posIdx)
        rhsStr = ' ';
    else
        parts = cell(numel(posIdx),1);
        for j = 1:numel(posIdx)
            parts{j} = sprintf('%g %s', full(S(posIdx(j),i)), mets{posIdx(j)});
        end
        rhsStr = strjoin(parts, ' + ');
    end
    rformulas{i} = sprintf('%s --> %s', lhsStr, rhsStr);
end

%% 4) Ensure tmp dir exists
if ~exist('tmp','dir'), mkdir('tmp'); end

%% 5) Configure EFMtool options
opts = CreateFluxModeOpts( ...
    'arithmetic',  'fractional', ...
    'compression', 'off',        ...
    'tmpdir',      'tmp',        ...
    'parallelize', true,         ...
    'timeout',     maxTimeSeconds ...
);

%% 6) Run EFMtool & catch timeout
fprintf('Starting EFM enumeration on %d reactions (timeout = %d s)…\n', nRxn, maxTimeSeconds);
t0 = tic;
try
    mnet = CalculateFluxModes(rformulas, rnames, opts);
    elapsed = toc(t0);
    nEFMs = size(mnet.efms, 2);
    fprintf('→ Returned after %.1f s, found %d EFMs\n', elapsed, nEFMs);
catch ME
    fprintf('EFMtool threw: %s\n', ME.message);
    partial = fullfile(opts.tmpdir, 'fluxModes.mat');
    if exist(partial,'file')
        tmp = load(partial,'efms');
        mnet.efms = tmp.efms;
        nEFMs     = size(mnet.efms,2);
        fprintf('→ Loaded partial %d EFMs from %s\n', nEFMs, partial);
    else
        rethrow(ME);
    end
end

% Save MATLAB results
save('efm_results.mat','mnet','opts','rnames');

%% 7) Export to CSV: reaction coefficients per mode
nEFMs = size(mnet.efms,2);
varNames = arrayfun(@(k) sprintf('Mode%d',k), 1:nEFMs,'Uni',false);
T = array2table(full(mnet.efms(:,1:nEFMs)), ...
       'VariableNames',varNames, 'RowNames',rnames);
writetable(T, '/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/ecoli/efmtool_coefficients.csv', 'WriteRowNames', true);
fprintf('Saved reaction‐by‐mode CSV with %d modes to efm_coefficients.csv\n', nEFMs);
%% 8) Print summary of EFMs in the same style as the MILP script

% Determine how many modes were actually found
nEFMs = size(mnet.efms, 2);

if restrictToSubsystem
    fprintf('\nFound %d EFMs in "%s":\n', nEFMs, subsystemName);
else
    fprintf('\nFound %d EFMs in full model:\n', nEFMs);
end

for k = 1:nEFMs
    % find which reactions participate
    supp = find(mnet.efms(:,k) ~= 0);
    reacList = rnames(supp);
    % print Mode number, support size, and reactions
    fprintf('  Mode %2d (|supp|=%d): %s\n', ...
            k, numel(supp), strjoin(reacList, ', '));
end
