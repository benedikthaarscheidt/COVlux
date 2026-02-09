%% 0) Paths & parameters
cobraModelFile = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/iML1515.mat';
outFile        = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/iML1515_pruned_permissive_biomassrm_loopless_thermo.mat';
changeCobraSolver('gurobi','LP');
changeCobraSolver('gurobi','MILP');
tol            = ((1e-5)-(1e-6));
coreBiomass    = 'BIOMASS_Ec_iML1515_core_75p37M';

keepBiomass      = false;
addDMfromBiomass = true;
removeBiomass    = true;

scaleBiomass   = false;
biomassNorm    = 'linf';
biomassTarget  = 1.0;

protectETC     = false;
protect_exact  = { ...
    'EX_o2_e','O2tex','O2tpp', ...
    'EX_co2_e','CO2tex', ...
    'NADH16','NADH18','NDH1','NDH2', ...
    'CYTBO3','CYTBD', ...
    'Q8','Q8H2', ...
    'ATPM', ...
    'BIOMASS_Ec_iML1515_core_75p37M', ...
};

usePermissiveMediumForPruning = true;

%% 1) Load the original reversible model
model_rev = readCbModel(cobraModelFile);
fprintf('1) Loaded reversible model: %d reactions\n\n', numel(model_rev.rxns));

rxns  = model_rev.rxns;
names = model_rev.rxnNames;
model_rev = changeObjective(model_rev, coreBiomass, 1);
bioIdx = find(strcmp(rxns, coreBiomass), 1);
if isempty(bioIdx)
    mask = ~cellfun('isempty', regexpi(rxns,  'biomass')) | ...
           ~cellfun('isempty', regexpi(names, 'biomass'));
    idx = find(mask);
    fprintf('Found %d biomass-like reactions (fallback search):\n', numel(idx));
    for k = idx(:).'
        fprintf('%4d  %-40s | %s\n', k, rxns{k}, names{k});
    end
    if isempty(idx)
        error('No biomass-like reaction found. Please set coreBiomass correctly.');
    end
    bioIdx = idx(1);
else
    fprintf('Found biomass by exact ID: %s (idx %d)\n', coreBiomass, bioIdx);
end

%% 2) Biomass handling: keep + scale (recommended), optional DM drains
bioCol     = model_rev.S(:,bioIdx);
substrates = find(bioCol < 0);

if addDMfromBiomass
    fprintf('2a) Adding DM_ drains for substrates of %s:\n', rxns{bioIdx});
    for i = 1:numel(substrates)
        r     = substrates(i);
        met   = model_rev.mets{r};
        safe  = regexprep(met,'[^A-Za-z0-9_]','_');
        newRxn = sprintf('DM_%s', safe);

        col = zeros(size(model_rev.S,1),1);
        col(r) = -1;

        model_rev.S(:,end+1)          = col;
        model_rev.rxns{end+1}         = newRxn;
        model_rev.lb(end+1,1)         = 0;
        model_rev.ub(end+1,1)         = 10000;
        model_rev.c(end+1,1)          = 0;
        model_rev.rules{end+1}        = '';
        model_rev.grRules{end+1}      = '';
        if isfield(model_rev,'rxnGeneMat')
            z = zeros(1, size(model_rev.rxnGeneMat,2));
            model_rev.rxnGeneMat(end+1,:) = z;
        end
        model_rev.rxnNames{end+1}     = newRxn;
        model_rev.subSystems{end+1}   = '';
        fprintf('   Added %s\n', newRxn);
    end
end

if keepBiomass && scaleBiomass
    col = model_rev.S(:,bioIdx);
    if strcmpi(biomassNorm,'linf')
        normCol = max(abs(col));
    elseif strcmpi(biomassNorm,'l1')
        normCol = sum(abs(col));
    else
        error('Unknown biomassNorm: %s', biomassNorm);
    end
    if normCol > 0
        s = normCol / biomassTarget;
        if abs(s - 1) > 1e-12
            model_rev.S(:,bioIdx) = model_rev.S(:,bioIdx) / s;
            if isfield(model_rev,'lb'), model_rev.lb(bioIdx) = model_rev.lb(bioIdx) * s; end
            if isfield(model_rev,'ub'), model_rev.ub(bioIdx) = model_rev.ub(bioIdx) * s; end
            if isfield(model_rev,'c'),  model_rev.c(bioIdx)  = model_rev.c(bioIdx)  / s; end
            fprintf('2b) Scaled biomass %s by 1/%.4g (norm=%s → %g)\n', ...
                    rxns{bioIdx}, s, biomassNorm, biomassTarget);
        else
            fprintf('2b) Biomass already at target norm (%s ≈ %g)\n', biomassNorm, biomassTarget);
        end
    end
end

%% 3) Prune blocked reactions from the reversible model

if usePermissiveMediumForPruning
    [model_rev, nEX, nUpt] = openMediumVerbose(model_rev);
    fprintf('Opened medium (REV): %d exchanges total, %d allow uptake (lb<0)\n\n', nEX, nUpt);
    dmMask = startsWith(string(model_rev.rxns), ["DM_","dm_"]);
    if any(dmMask)
        model_rev.lb(dmMask) = max(model_rev.lb(dmMask), 0);
        if isfield(model_rev,'rev'), model_rev.rev(dmMask) = 0; end
        idx = find(dmMask);
        fprintf('Re-closed DM_ drains: %d reactions. Examples:\n', nnz(dmMask));
        disp(model_rev.rxns(idx(1:min(10,end))));
    end
    sinkMask = startsWith(string(model_rev.rxns), ["SINK_","sink_","SK_","sk_"]);
    model_rev.lb(sinkMask) = max(model_rev.lb(sinkMask), 0);
    if isfield(model_rev,'rev'), model_rev.rev(sinkMask) = 0; end
end

iter = 0;
while true
    iter = iter + 1;
    fprintf('3.%d) Checking blocked rxns in reversible model (fluxVariability)...\n', iter);
    [minFlux, maxFlux] = fluxVariability(model_rev, 1, 'max', model_rev.rxns, 0, 0);
    blocked = (maxFlux < tol) & (minFlux > -tol);
    toRemove = find(blocked);
    if isempty(toRemove)
        fprintf('3.%d) No more blocked reactions (rev).\n\n', iter);
        break;
    end
    if protectETC
        keepMask = false(size(toRemove));
        for ii = 1:numel(toRemove)
            rname = model_rev.rxns{toRemove(ii)};
            if any(strcmp(rname, protect_exact))
                keepMask(ii) = true;
            end
        end
        if any(keepMask)
            kept = toRemove(keepMask);
            fprintf('   Keeping %d protected rxns (rev):\n', numel(kept));
            for k = kept(:).'
                fprintf('     * %s\n', model_rev.rxns{k});
            end
        end
        toRemove = toRemove(~keepMask);
    end
    if isempty(toRemove)
        fprintf('   Nothing left to remove this round (rev).\n\n');
        break;
    end
    toRemove = sort(toRemove,'descend');
    fprintf('   Removing %d blocked rxns (rev):\n', numel(toRemove));
    for k = toRemove.'
        fprintf('     - %s\n', model_rev.rxns{k});
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
model_ir = changeObjective(model_ir, coreBiomass, 1);

%% 5) Prune blocked reactions from the irreversible model using MILP with f/b ban and boundary activation

iter = 0;
while true
    iter = iter + 1;
    fprintf('5.%d) Checking blocked rxns in irreversible model (MILP-FVA)...\n', iter);
    maxFlux_ir = fvaMaxMILP_pairBan_and_boundary(model_ir, tol);
    blocked = maxFlux_ir < tol;
    toRemove = find(blocked);
    if isempty(toRemove)
        fprintf('5.%d) No more blocked reactions (irrev).\n\n', iter);
        break;
    end
    if protectETC
        keepMask = false(size(toRemove));
        for ii = 1:numel(toRemove)
            rname = model_ir.rxns{toRemove(ii)};
            if any(strcmp(rname, protect_exact))
                keepMask(ii) = true;
            end
        end
        if any(keepMask)
            kept = toRemove(keepMask);
            fprintf('   Keeping %d protected rxns (irrev):\n', numel(kept));
            for k = kept(:).'
                fprintf('     * %s\n', model_ir.rxns{k});
            end
        end
        toRemove = toRemove(~keepMask);
    end
    if isempty(toRemove)
        fprintf('   Nothing left to remove this round (irrev).\n\n');
        break;
    end
    toRemove = sort(toRemove,'descend');
    fprintf('   Removing %d blocked rxns (irrev):\n', numel(toRemove));
    for k = toRemove.'
        fprintf('     - %s\n', model_ir.rxns{k});
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

rxns=model_ir.rxns;
bioIdx = find(strcmp(rxns, coreBiomass), 1);
if removeBiomass && ~keepBiomass
    fprintf('2c) Removing core biomass %s (idx %d)\n\n', rxns{bioIdx}, bioIdx);
    model_ir.S(:,        bioIdx)        = [];
    model_ir.rxns(      bioIdx)        = [];
    model_ir.lb(        bioIdx)        = [];
    model_ir.ub(        bioIdx)        = [];
    if isfield(model_ir,'c'),          model_ir.c(bioIdx)          = []; end
    if isfield(model_ir,'rules'),      model_ir.rules(bioIdx)      = []; end
    if isfield(model_ir,'grRules'),    model_ir.grRules(bioIdx)    = []; end
    if isfield(model_ir,'rxnGeneMat'), model_ir.rxnGeneMat(bioIdx,:) = []; end
    if isfield(model_ir,'rxnNames'),   model_ir.rxnNames(bioIdx)   = []; end
    if isfield(model_ir,'subSystems'), model_ir.subSystems(bioIdx) = []; end
elseif removeBiomass && keepBiomass
    warning('removeBiomass=true but keepBiomass=true; ignoring removal and keeping Biomass.');
end

model_ir = removeOrphanMets_verbose(model_ir);

pruned_ir = model_ir;

save(outFile, 'pruned_ir');
fprintf('6) Saved pruned & split model:\n   %s\n', outFile);
fprintf('   Total reactions: %d\n', numel(pruned_ir.rxns));

%%

function [model, nEX, nUpt] = openMediumVerbose(model)
    ex = findExcRxns(model);
    nEX = nnz(ex);
    model.lb(ex) = min(model.lb(ex), -1000);
    model.ub(ex) = max(model.ub(ex),  1000);
    nUpt = nnz(model.lb(ex) < 0);
    closed = find(ex & ~(model.lb < 0));
    if ~isempty(closed)
        fprintf('Examples of exchanges with no uptake (first 10):\n');
        disp(model.rxns(closed(1:min(10,end))));
    end
end

function [model, masks] = openMedium_REV_viaCOBRA(model, EXcap, DMcap)
    if nargin<2 || isempty(EXcap), EXcap = 1000; end
    if nargin<3 || isempty(DMcap), DMcap = 1000; end
    [~, ExchRxnBool, DMRxnBool, SinkRxnBool] = findSExRxnInd(model);
    exOnly = ExchRxnBool & ~DMRxnBool & ~SinkRxnBool;
    if any(DMRxnBool)
        model.lb(DMRxnBool) = max(model.lb(DMRxnBool), 0);
        model.ub(DMRxnBool) = max(model.ub(DMRxnBool), DMcap);
        if isfield(model,'rev'), model.rev(DMRxnBool) = 0; end
    end
    if any(SinkRxnBool)
        model.lb(SinkRxnBool) = max(model.lb(SinkRxnBool), 0);
        model.ub(SinkRxnBool) = max(model.ub(SinkRxnBool), DMcap);
        if isfield(model,'rev'), model.rev(SinkRxnBool) = 0; end
    end
    if any(exOnly)
        model.lb(exOnly) = min(model.lb(exOnly), -EXcap);
        model.ub(exOnly) = max(model.ub(exOnly),  EXcap);
    end
    masks.exOnly      = exOnly;
    masks.dmBool      = DMRxnBool;
    masks.sinkBool    = SinkRxnBool;
    masks.nEX_opened  = nnz(exOnly);
    masks.nUptakeOpen = nnz(model.lb(exOnly) < 0);
    fprintf('openMedium_REV_viaCOBRA: opened %d EX (uptake on for %d)\n', ...
            masks.nEX_opened, masks.nUptakeOpen);
    if masks.nEX_opened>0
        exIdx = find(exOnly);
        disp(model.rxns(exIdx(1:min(10,end))));
    end
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

function maxFlux = fvaMaxMILP_pairBan_and_boundary(model, epsBoundary, alphaGrowth)
    if nargin<2 || isempty(epsBoundary), epsBoundary = 1e-6; end
    if nargin<3 || isempty(alphaGrowth), alphaGrowth = 0.01; end
    n = numel(model.rxns);
    S = model.S; lb = model.lb(:); ub = model.ub(:);
    exBool = findExcRxns(model);
    rxs = string(model.rxns);
    dmMask = startsWith(rxs,["DM_","dm_"]);
    sinkMask = startsWith(rxs,["SINK_","sink_","SK_","sk_"]);
    exOnly = exBool & ~dmMask & ~sinkMask;
    baseNames = regexprep(model.rxns,'(_f|_b)$','');
    isF = ~cellfun('isempty', regexp(model.rxns,'_f$'));
    isB = ~cellfun('isempty', regexp(model.rxns,'_b$'));
    [~,~,ib] = unique(baseNames);
    pairs = [];
    ug = unique(ib(:)).';
    for g = ug
        idx = find(ib==g);
        if numel(idx)>=2
            f = idx(isF(idx)); b = idx(isB(idx));
            if ~isempty(f) && ~isempty(b), pairs = [pairs; [f(1) b(1)]]; end
        end
    end
    mPairs = size(pairs,1);
    totalVars = n + 2*mPairs;
    A = sparse(0,totalVars); bvec = []; csense = '';
    Ae = [S sparse(size(S,1),2*mPairs)]; be = zeros(size(S,1),1); csE = repmat('E',size(S,1),1);
    A = [A; Ae]; bvec = [bvec; be]; csense = [csense; csE];
    rows = 3*mPairs; Afb = spalloc(rows,totalVars,rows*3); bf = zeros(rows,1); csfb = repmat('L',rows,1); r = 0;
    for p = 1:mPairs
        f = pairs(p,1); bb = pairs(p,2);
        Mf = max(0,ub(f)); Mb = max(0,ub(bb));
        r = r + 1; Afb(r,f)  = 1; Afb(r,n + 2*p - 1) = -Mf;
        r = r + 1; Afb(r,bb) = 1; Afb(r,n + 2*p    ) = -Mb;
        r = r + 1; Afb(r,n + 2*p - 1) = 1; Afb(r,n + 2*p) = 1; bf(r) = 1;
    end
    A = [A; Afb]; bvec = [bvec; bf]; csense = [csense; csfb];
    A = [A; sparse(1,totalVars)]; A(end,1:n) = double(exOnly(:))'; bvec = [bvec; max(epsBoundary,1e-9)]; csense = [csense; 'G'];
    bioIdx = find(model.c~=0,1);
    if ~isempty(bioIdx) && alphaGrowth>0
        solFBA = optimizeCbModel(model,'max');
        if ~isempty(solFBA.f) && solFBA.f>0
            row = sparse(1,totalVars); row(1,bioIdx) = 1;
            A = [A; row]; bvec = [bvec; alphaGrowth*solFBA.f]; csense = [csense; 'G'];
        end
    end
    lbMILP = [lb; zeros(2*mPairs,1)];
    ubMILP = [ub; ones(2*mPairs,1)];
    vartype = [repmat('C',1,n) repmat('B',1,2*mPairs)];
    baseProb.A = A; baseProb.b = bvec; baseProb.csense = csense; baseProb.lb = lbMILP; baseProb.ub = ubMILP; baseProb.vartype = vartype; baseProb.osense = -1;
    maxFlux = zeros(n,1);
    for i = 1:n
        c = zeros(totalVars,1); c(i) = 1; baseProb.c = c;
        sol = solveCobraMILP(baseProb);
        xi = 0;
        if isstruct(sol)
            if isfield(sol,'full') && ~isempty(sol.full), xi = sol.full(i);
            elseif isfield(sol,'x') && ~isempty(sol.x),   xi = sol.x(i);
            elseif isfield(sol,'obj') && ~isempty(sol.obj), xi = sol.obj;
            end
        end
        maxFlux(i) = xi;
    end
end
