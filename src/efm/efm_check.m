%% verifyEFMs_fullPipeline.m
% 1) Load & optionally restrict to subsystem
% 2) Convert to irreversible
% 3) Load your CSV of “Mode#” flux‐vectors
% 4) For each mode:
%       – lb(i)=eps_flux, ub(i)=original_ub  if i ∈ support
%       – lb=ub=0                             if i ∉ support
%       – solve S·v=0, v≥0 with linprog (zero obj)
%       – if feasible, test minimality by dropping each i∈support

%% ——— USER SETUP ———————————————————————————————
cobraModelFile      = 'Recon3D_301.mat';
restrictToSubsystem = false;    % ← match what you did when enumerating
subsystemName       = 'Glycolysis/Gluconeogenesis';
csvFile             = 'efm_fluxes_output.csv';
eps_flux            = 1e-3;    % same epsilon you used in the MILP
maxFluxUB           = 1e3;    % cap for any ub to keep LP bounded
tol                 = 1e-6;   % numerical tolerance
%% ——————————————————————————————————————————————

%% 1) Load & slice
orig = readCbModel(cobraModelFile);
if restrictToSubsystem
    ss = orig.subSystems;
    % normalize to a cellstr of names
    if ischar(ss)
        ss = cellstr(ss);
    elseif iscell(ss) && all(cellfun(@iscell, ss))
        ss = cellfun(@(c) strjoin(c,'; '), ss, 'Uni',false);
    end
    mask = contains(ss, subsystemName, 'IgnoreCase',true);
    if ~any(mask)
        error('No reactions in subsystem "%s".', subsystemName);
    end
    fprintf('Restricting to %d reactions in "%s".\n', nnz(mask), subsystemName);
    model.S    = orig.S(:, mask);
    model.rxns = orig.rxns(mask);
    model.lb   = orig.lb(mask);
    model.ub   = orig.ub(mask);
    if isfield(orig,'c'), model.c = orig.c(mask); end
else
    fprintf('Using full model with %d reactions.\n', numel(orig.rxns));
    model = orig;
end

%% 2) Irreversible conversion
model_ir = convertToIrreversible(model);
S        = model_ir.S;
nRxns    = numel(model_ir.rxns);

%% 3) Read your CSV of modes
T          = readtable(csvFile, 'ReadRowNames',true,'VariableNamingRule','preserve');
modeNames  = T.Properties.VariableNames;
rxnCSV     = T.Properties.RowNames;
nModes     = numel(modeNames);

%% 4) Prepare LP options & base bounds
opts = optimoptions('linprog','Display','none');
lb_base = zeros(nRxns,1);
ub_base = model_ir.ub;
ub_base(isinf(ub_base)) = maxFluxUB;

results = struct();

for k = 1:nModes
    mode   = modeNames{k};
    vCSV   = T.(mode);
    
    % 4a) Determine support exactly as y==0 in your MILP
    suppMask = abs(vCSV) >= eps_flux;
    suppRxns  = rxnCSV(suppMask);
    fprintf('\nMode %s: %d reactions in support\n', mode, sum(suppMask));
    
    % map onto model_ir indices
    isSupp = ismember(model_ir.rxns, suppRxns);
    
    % 4b) Build bounds:  on-support [eps_flux, orig_ub], off-support = 0
    lb = lb_base;  ub = ub_base;
    lb(isSupp) = eps_flux;
    ub(~isSupp) = 0;
    
    % 5) Feasibility via S·v=0
    [~,~,exitflag] = linprog(zeros(nRxns,1),[],[], S, zeros(size(S,1),1), lb, ub, opts);
    feas = (exitflag == 1);
    fprintf('  Feasibility: %s\n', tern(feas,'PASS','FAIL'));
    
    % 6) Minimality: drop each reaction and re‐solve
    minimal = true;
    if feas
        for j = find(isSupp)'
            lb2 = lb; ub2 = ub;
            lb2(j) = 0; ub2(j) = 0;  % knock out reaction j
            
            [~,~,ef2] = linprog(zeros(nRxns,1),[],[], S, zeros(size(S,1),1), lb2, ub2, opts);
            if ef2 == 1
                fprintf('    Removing %s still feasible → NOT minimal\n', model_ir.rxns{j});
                minimal = false;
                break;
            end
        end
    else
        minimal = false;
    end
    
    fprintf('  Minimality: %s\n', tern(minimal,'PASS','FAIL'));
    results.(mode) = struct('feasible',feas,'minimal',minimal);
end

%% 7) Summary
nTrue = sum(structfun(@(r) r.feasible && r.minimal, results));
fprintf('\n%d / %d supports are true minimal EFMs.\n', nTrue, nModes);

function out = tern(cond,a,b)
    if cond, out = a; else out = b; end
end
