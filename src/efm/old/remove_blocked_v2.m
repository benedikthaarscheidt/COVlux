%% remove_blocked_and_split_biomass_v11.m
% Two‐stage FVA‐based pruning (reversible → irreversible), manual reaction removal,
% split out the core biomass sinks up front, removing entries from all relevant fields.

%% 0) Paths & parameters
cobraModelFile = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/iML1515.mat';
outFile        = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/iML1515_final.mat';
tol            = 1e-5;
coreBiomass    = "BIOMASS_Ec_iML1515_core_75p37M";
%'Biomass Objective Function with GAM';%'BIOMASS_Ec_iML1515_core_75p37M';

%% 1) Load the original reversible model
model_rev = readCbModel(cobraModelFile);
fprintf('1) Loaded reversible model: %d reactions\n\n', numel(model_rev.rxns));

rxns  = model_rev.rxns;
names = model_rev.rxnNames;

% Case-insensitive match for the word "biomass" in either IDs or names
mask = ~cellfun('isempty', regexpi(rxns,  'biomass')) | ...
       ~cellfun('isempty', regexpi(names, 'biomass'));

idx = find(mask);

% Print nicely
fprintf('Found %d biomass-like reactions:\n', numel(idx));
for k = idx(:).'
    fprintf('%4d  %-40s | %s\n', k, rxns{k}, names{k});
end
bioIdx = idx(1);
%% 2) Split & remove the core biomass reaction

if isempty(bioIdx)
    error('Core biomass reaction %s not found.', coreBiomass);
end

bioCol     = model_rev.S(:,bioIdx);
substrates = find(bioCol < 0);

fprintf('2) Splitting core biomass (%s) into DM_ drains:\n', coreBiomass);
for i = 1:numel(substrates)
    r     = substrates(i);
    coeff = bioCol(r);                     % negative coefficient
    met   = model_rev.mets{r};
    safe  = regexprep(met,'[^A-Za-z0-9_]','_');
    newRxn = sprintf('DM_%s', safe);

    col = zeros(size(model_rev.S,1),1);
    col(r) = -1;                       

    % append new sink reaction
    model_rev.S(:,end+1)         = col;
    model_rev.rxns{end+1}        = newRxn;
    model_rev.lb(end+1,1)        = 0;
    model_rev.ub(end+1,1)        = Inf;
    model_rev.c(end+1,1)         = 0;

    % copy all gene-associated fields from the biomass reaction
    model_rev.rules{end+1}       = model_rev.rules{bioIdx};
    model_rev.grRules{end+1}     = model_rev.grRules{bioIdx};
    model_rev.rxnGeneMat(end+1,:) = model_rev.rxnGeneMat(bioIdx,:);

   
    model_rev.rxnNames{end+1}    = newRxn;
    model_rev.subSystems{end+1}  = '';

    fprintf('   Added %s (stoich %.4g)\n', newRxn, coeff);
end

fprintf('   Removing core biomass %s (idx %d)\n\n', coreBiomass, bioIdx);
% manual removal from all fields
model_rev.S(:,        bioIdx)        = [];
model_rev.rxns(      bioIdx)        = [];
model_rev.lb(        bioIdx)        = [];
model_rev.ub(        bioIdx)        = [];
if isfield(model_rev,'c'),          model_rev.c(bioIdx)          = []; end
if isfield(model_rev,'rules'),      model_rev.rules(bioIdx)      = []; end
if isfield(model_rev,'grRules'),    model_rev.grRules(bioIdx)    = []; end
if isfield(model_rev,'rxnGeneMat'), model_rev.rxnGeneMat(bioIdx,:) = []; end
if isfield(model_rev,'rxnNames'),   model_rev.rxnNames(bioIdx)   = []; end
if isfield(model_rev,'subSystems'), model_rev.subSystems(bioIdx) = []; end

%% stoichiometric adjustment to obtain a full integer model.S
%tol = 1e-8;
%S = model_rev.S;
%[nMets, nRxns] = size(S);
%scaled_rxns = [];
%
%for j = 1:nRxns
%    sj = S(:,j);
%    nonzero = find(abs(sj) > tol);
%    if isempty(nonzero), continue; end
%    
%    vals = sj(nonzero);
%    int_diff = abs(vals - round(vals));
%    if any(int_diff > tol)
%        % Find the *smallest* non-integer magnitude value
%        frac_vals = vals(int_diff > tol);
%        min_val = min(abs(frac_vals));
%        factor = 1 / min_val;
%        % Round factor to nearest integer to avoid floating point artifacts
%        factor = round(factor);
%        S(:,j) = S(:,j) * factor;
%        if isfield(model_rev, 'lb'), model_rev.lb(j) = model_rev.lb(j) / factor; end
%        if isfield(model_rev, 'ub'), model_rev.ub(j) = model_rev.ub(j) / factor; end
%        scaled_rxns(end+1) = j;
%    end
%end
%
%model_rev.S = S;
%
%fprintf('\nScaled %d reactions to integer stoichiometry:\n', numel(scaled_rxns));
%for k = 1:numel(scaled_rxns)
%    idx = scaled_rxns(k);
%    fprintf('%4d: %s\n', idx, model_rev.rxns{idx});
%end
%% 3) Prune blocked reactions from the reversible model. We do this in the reversible model as when we only remove in the irreversible model 
% a reversible reaction which was split into _f and _b and blocked priorly can suddenly have flux as both _f and _b can be active at the same time (
% which is wrong and not removing blocked reactions properly)

iter = 0;
while true
    iter = iter + 1;
    fprintf('3.%d) Checking blocked rxns in reversible model (fluxVariability)...\n', iter);

    [minFlux, maxFlux] = fluxVariability(model_rev, tol, ...
                                         'rxnNameList', model_rev.rxns);
    blocked = (maxFlux < tol) & (minFlux > -tol);

    toRemove = find(blocked);
    if isempty(toRemove)
        fprintf('3.%d) No more blocked reactions (rev).\n\n', iter);
        break;
    end

    toRemove = sort(toRemove,'descend');
    fprintf('   Removing %d blocked rxns:\n', numel(toRemove));
    for k = toRemove'
        fprintf('     - %s\n', model_rev.rxns{k});
        % manual removal from all fields
        model_rev.S(:,        k)        = [];
        model_rev.rxns(      k)        = [];
        model_rev.lb(        k)        = [];
        model_rev.ub(        k)        = [];
        if isfield(model_rev,'c'),          model_rev.c(k)          = []; end
        if isfield(model_rev,'rules'),      model_rev.rules(k)      = []; end
        if isfield(model_rev,'grRules'),    model_rev.grRules(k)    = []; end
        if isfield(model_rev,'rxnGeneMat'), model_rev.rxnGeneMat(k,:) = []; end
        if isfield(model_rev,'rxnNames'),   model_rev.rxnNames(k)   = []; end
        if isfield(model_rev,'subSystems'), model_rev.subSystems(k) = []; end
    end
    fprintf('\n');
end

%% 4) Convert the pruned reversible model to irreversible
model_ir = convertToIrreversible(model_rev);
fprintf('4) Converted to irreversible: %d reactions\n\n', numel(model_ir.rxns));

%% 5) Prune blocked reactions from the irreversible model
iter = 0;
while true
    iter = iter + 1;
    fprintf('5.%d) Checking blocked rxns in irreversible model (fluxVariability)...\n', iter);

    [~, maxFlux_ir] = fluxVariability(model_ir, tol, ...
                                      'rxnNameList', model_ir.rxns);
    blocked = maxFlux_ir < tol;

    toRemove = find(blocked);
    if isempty(toRemove)
        fprintf('5.%d) No more blocked reactions (irrev).\n\n', iter);
        break;
    end

    toRemove = sort(toRemove,'descend');
    fprintf('   Removing %d blocked rxns:\n', numel(toRemove));
    for k = toRemove'
        fprintf('     - %s\n', model_ir.rxns{k});
        % manual removal from all fields
        model_ir.S(:,        k)        = [];
        model_ir.rxns(      k)        = [];
        model_ir.lb(        k)        = [];
        model_ir.ub(        k)        = [];
        if isfield(model_ir,'c'),          model_ir.c(k)          = []; end
        if isfield(model_ir,'rules'),      model_ir.rules(k)      = []; end
        if isfield(model_ir,'grRules'),    model_ir.grRules(k)    = []; end
        if isfield(model_ir,'rxnGeneMat'), model_ir.rxnGeneMat(k,:) = []; end
        if isfield(model_ir,'rxnNames'),   model_ir.rxnNames(k)   = []; end
        if isfield(model_ir,'subSystems'), model_ir.subSystems(k) = []; end
    end
    fprintf('\n');
end


%% remove orphan metabolites

model_ir = removeOrphanMets(model_ir);

%% 6) Save the final model
pruned_ir = model_ir;
%%% this is for making non integer entries in S integer by scaling the reactions coefficients. This didnt help in resolving the isses and can therefor be ignored.
%S = model_ir.S;
%rxns = model_ir.rxns;
%

%nonint_rxns = [];
%tol = 1e-8;  % tolerance for floating point errors
%
%for j = 1:size(S,2)
%    sj = S(:,j);
%    if any(abs(sj - round(sj)) > tol & sj ~= 0)
%        nonint_rxns(end+1) = j;
%    end
%end
%
%fprintf('\nReactions with non-integer stoichiometry (%d):\n', numel(nonint_rxns));
%for k = 1:numel(nonint_rxns)
%    idx = nonint_rxns(k);
%    fprintf('%4d: %s\n', idx, rxns{idx});
%end

%% From here on it is only for checking the model to see if there are any disconnected components. 

S = model_ir.S;
rxns = model_ir.rxns;
mets = model_ir.mets;


[m, n] = size(S);
A = [ sparse(m, m),        abs(S) > 0;
      (abs(S) > 0)',   sparse(n, n) ];


G = graph(A | A');
component = conncomp(G);


met_components = component(1:m);
rxn_components = component(m+1:end);


num_components = max(component);
counts = histcounts(component, 1:num_components+1);

fprintf('\nNetwork connectivity analysis:\n');
for k = 1:num_components
    mets_in = find(met_components == k);
    rxns_in = find(rxn_components == k);
    fprintf('  Component %d: %d metabolites, %d reactions\n', ...
        k, numel(mets_in), numel(rxns_in));
end

% List disconnected reactions (not in the largest component)
[~,biggest] = max(counts);
main_component = biggest;
disconnected_rxns = find(rxn_components ~= main_component);
disconnected_mets = find(met_components ~= main_component);

fprintf('\nDisconnected reactions (%d):\n', numel(disconnected_rxns));
for i = 1:numel(disconnected_rxns)
    idx = disconnected_rxns(i);
    fprintf('%4d: %s\n', idx, rxns{idx});
end
fprintf('\nDisconnected metabolites (%d):\n', numel(disconnected_mets));
for i = 1:numel(disconnected_mets)
    idx = disconnected_mets(i);
    fprintf('%4d: %s\n', idx, mets{idx});
end
%%


save(outFile, 'pruned_ir');
fprintf('6) Saved pruned & split model:\n   %s\n', outFile);
fprintf('   Total reactions: %d\n', numel(pruned_ir.rxns));


%%

function model = removeOrphanMets(model)
    deg = full(sum(abs(model.S)>0,2));
    toRemove = find(deg==0);
    if isempty(toRemove)
        return
    end
    fprintf('Removing %d orphan metabolites:\n', numel(toRemove));
    for i = 1:numel(toRemove)
        fprintf('  - %s\n', model.mets{toRemove(i)});
    end
    model.S(toRemove,:) = [];
    if isfield(model,'mets'),            model.mets(toRemove)            = []; end
    if isfield(model,'metNames'),         model.metNames(toRemove)         = []; end
    if isfield(model,'metFormulas'),      model.metFormulas(toRemove)      = []; end
    if isfield(model,'metCharge'),        model.metCharge(toRemove)        = []; end
    if isfield(model,'metComps'),         model.metComps(toRemove)         = []; end
    if isfield(model,'metCompartment'),   model.metCompartment(toRemove)   = []; end
    if isfield(model,'b'),                model.b(toRemove)                = []; end
    if isfield(model,'csense')
        if ischar(model.csense)
            model.csense(toRemove) = [];
        else
            model.csense(toRemove) = model.csense(setdiff(1:numel(model.csense),toRemove));
        end
    end
    if isfield(model,'metKEGGID'),        model.metKEGGID(toRemove)        = []; end
    if isfield(model,'metHMDBID'),        model.metHMDBID(toRemove)        = []; end
    if isfield(model,'metChEBIID'),       model.metChEBIID(toRemove)       = []; end
    if isfield(model,'metPubChemID'),     model.metPubChemID(toRemove)     = []; end
    if isfield(model,'metInChIString'),   model.metInChIString(toRemove)   = []; end
    if isfield(model,'metNotes'),         model.metNotes(toRemove)         = []; end
end

function model = removeOrphanMets_verbose(model)
    deg = full(sum(abs(model.S)>0,2));
    toRemove = find(deg==0);
    if isempty(toRemove)
        fprintf('No orphan metabolites to remove.\n');
        return
    end
    fprintf('Removing %d orphan metabolites:\n', numel(toRemove));
    for i = 1:numel(toRemove)
        fprintf('  - %s\n', model.mets{toRemove(i)});
    end
    model.S(toRemove,:) = [];
    if isfield(model,'mets'),            model.mets(toRemove)            = []; end
    if isfield(model,'metNames'),        model.metNames(toRemove)         = []; end
    if isfield(model,'metFormulas'),     model.metFormulas(toRemove)      = []; end
    if isfield(model,'metCharge'),       model.metCharge(toRemove)        = []; end
    if isfield(model,'metComps'),        model.metComps(toRemove)         = []; end
    if isfield(model,'metCompartment'),  model.metCompartment(toRemove)   = []; end
    if isfield(model,'b'),               model.b(toRemove)                = []; end
    if isfield(model,'csense')
        if ischar(model.csense)
            model.csense(toRemove) = [];
        else
            model.csense(toRemove) = model.csense(setdiff(1:numel(model.csense),toRemove));
        end
    end
    if isfield(model,'metKEGGID'),       model.metKEGGID(toRemove)        = []; end
    if isfield(model,'metHMDBID'),       model.metHMDBID(toRemove)        = []; end
    if isfield(model,'metChEBIID'),      model.metChEBIID(toRemove)       = []; end
    if isfield(model,'metPubChemID'),    model.metPubChemID(toRemove)     = []; end
    if isfield(model,'metInChIString'),  model.metInChIString(toRemove)   = []; end
    if isfield(model,'metNotes'),        model.metNotes(toRemove)         = []; end
end

%%



[model_rev, nEX, nUpt] = openMediumVerbose(model_rev);
fprintf('Opened medium for REV pruning: %d exchanges total, %d allow uptake (lb<0)\n\n', nEX, nUpt);