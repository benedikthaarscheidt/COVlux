% =========================================================================
% MRAS PIPELINE: CLUSTERED DATA -> REACTION ACTIVITY & COVARIANCE
% Handles: Gene Symbols, Operon Splitting, and Locus Tags
% =========================================================================

% ---- 1. SETUP & CONFIG --------------------------------------------------
currentScriptPath = fileparts(mfilename('fullpath'));
projectRoot       = fileparts(fileparts(currentScriptPath)); % src / mras -> src -> root
configFile = fullfile(projectRoot, 'config', 'config.json');


if ~isfile(configFile)
    error('Config file not found: %s', configFile);
end

config = jsondecode(fileread(configFile));

% ---- 2. DYNAMIC PATH RESOLUTION -----------------------------------------
% Model path from config
modelPath = fullfile(projectRoot, config.paths.models_dir, config.model.model_file);

% Determine Input Directory (dataDir) based on config
clustering_base = fullfile(projectRoot, config.paths.results_dir, 'Clustering');
run_name        = config.params.input_clustering_folder;
second_moment_mat = config.params.second_moment_mat;
if isfield(config.params, 'use_conditions') && config.params.use_conditions
    dataDir = fullfile(clustering_base, run_name, 'data_files', 'grouped_by_condition');
    if ~exist(dataDir, 'dir')
        fprintf('Warning: Grouped directory not found, falling back to standard data_files: %s\n', dataDir);
        dataDir = fullfile(clustering_base, run_name, 'data_files');
    end
else
    dataDir = fullfile(clustering_base, run_name, 'data_files');
end

outDir = fullfile(dataDir, 'MRAS_outputs');
fprintf('MRAS Input Directory: %s\n', dataDir);
fprintf('MRAS Output Directory: %s\n', outDir);

% ---- 3. CONFIGURATION & CONTROLS ----------------------------------------
% TRUE  = Augment map with Gene Names/Symbols (dnaK -> b0014) AND split operons (A-B -> A, B)
% FALSE = Strict matching only (Data header MUST match model.genes exactly)
use_smart_mapping = config.params.use_smart_mapping;

% EXPRESSION SCALE CONTROL
input_scale = 'log2';
log_offset  = 1.0;

% -------------------------------------------------------------------------
if ~exist(outDir, 'dir'), mkdir(outDir); end

% ---- 4. LOAD MODEL (ONCE)
fprintf('Loading model...\n');
S = load(modelPath); 
fn = fieldnames(S); 
model = S.(fn{1});
genes   = string(model.genes(:)); 
rxns    = string(model.rxns(:)); 
grRules = string(model.grRules(:));

% ---- 5. BUILD GENE MAP (SMART vs SIMPLE)
fprintf('Building Gene Index Map (Smart Mode: %d)...\n', use_smart_mapping);
geneIndexMap = containers.Map();
% A. Always map standard IDs (b0001 -> 1)
for i = 1:numel(genes)
    geneIndexMap(lower(genes(i))) = i;
end
% B. Conditional: Map Names/Symbols (dnaK -> 1)
if use_smart_mapping
    if isfield(model, 'geneNames')
        fprintf('  Augmenting map with model.geneNames...\n');
        gNames = string(model.geneNames(:));
        for i = 1:numel(gNames)
            key = lower(gNames(i));
            if strlength(key) > 0 && ~isKey(geneIndexMap, key)
                geneIndexMap(key) = i;
            end
        end
    end
    if isfield(model, 'geneSymbols')
        fprintf('  Augmenting map with model.geneSymbols...\n');
        gSymb = string(model.geneSymbols(:));
        for i = 1:numel(gSymb)
            key = lower(gSymb(i));
            if strlength(key) > 0 && ~isKey(geneIndexMap, key)
                geneIndexMap(key) = i;
            end
        end
    end
end
fprintf('  Gene Map ready with %d unique keys.\n', geneIndexMap.Count);

% ---- 6. COMPUTE GENE PROMISCUITY (ONCE)
fprintf('Computing gene promiscuity...\n');
if isfield(model,'rxnGeneMat')
    gcount = sum(model.rxnGeneMat~=0, 1); 
    gcount = gcount(:)'; 
    if numel(gcount) ~= numel(genes)
        gcount = count_gene_promiscuity_text(cellstr(grRules), cellstr(genes)); 
    end
else
    gcount = count_gene_promiscuity_text(cellstr(grRules), cellstr(genes));
end
gcount = max(gcount, 1);

% ---- 7. PARSE GPR RULES (ONCE)
fprintf('Parsing GPR rules...\n');
postfixAll = cell(size(grRules));
for j = 1:numel(grRules)
    r = grRules(j);
    if strlength(strtrim(r)) == 0
        postfixAll{j} = {};
    else
        postfixAll{j} = gpr_to_postfix(char(r), geneIndexMap);
    end
end

% ---- 8. FILE DISCOVERY
files = dir(fullfile(dataDir, 'cluster_*')); 
files = files(~[files.isdir]);
mask = endsWith({files.name}, '.csv', 'IgnoreCase', true); 
files = files(mask);
[~, uIdx] = unique(lower(fullfile({files.folder}, {files.name}))); 
files = files(uIdx);
fprintf('Found %d cluster files\n', numel(files));

% ---- 9. MAIN PROCESSING LOOP
for f = 1:numel(files)
    tic;
    infile = fullfile(files(f).folder, files(f).name); 
    fprintf('\n=== Processing %s ===\n', files(f).name);
    
    % A. Read Table
    T = readtable(infile, 'VariableNamingRule', 'preserve', 'TextType', 'string');
    numVars = varfun(@isnumeric, T, 'OutputFormat', 'uniform');
    Xraw = table2array(T(:, numVars)); 
    usedNames = T.Properties.VariableNames(numVars);
    
    fprintf('Input data: %d cells x %d genes\n', size(Xraw, 1), size(Xraw, 2));
    
    % B. Scale Conversion
    switch lower(input_scale)
        case 'linear', Xlin = Xraw;
        case 'log2',   Xlin = 2.^Xraw - log_offset; Xlin(Xlin<0) = 0;
        case 'ln',     Xlin = exp(Xraw) - log_offset; Xlin(Xlin<0) = 0;
    end

    row_var = var(Xraw, 0, 2, 'omitnan');
    if mean(row_var) < 1e-6
        warning('!!! INPUT DATA IS FLAT !!! The CSV file itself has identical values per cell.');
    else
        fprintf('[Pass] Input CSV has healthy variance (Mean Var: %.4f)\n', mean(row_var));
    end
    
    % C. Collapse Duplicate Columns
    L = lower(string(usedNames)); 
    [uniqNames, ~, ic] = unique(L, 'stable'); 
    fprintf('After duplicate collapse: %d unique genes\n', numel(uniqNames));
    
    Xuniq = zeros(size(Xlin, 1), numel(uniqNames));
    for k = 1:numel(uniqNames)
        cols = (ic == k); 
        if sum(cols) > 1
            Xuniq(:, k) = mean(Xlin(:, cols), 2, 'omitnan');
        else
            Xuniq(:, k) = Xlin(:, cols);
        end
    end
    
   
    % =========================================================================
    % D. MAP DATA TO MODEL (Diagnostic Mode)
    % =========================================================================
    Gmat = zeros(size(Xuniq, 1), numel(genes)); 
    
    % Initialize BaseName and Accumulators
    [~, baseName, ~] = fileparts(files(f).name);
    Gmat_acc   = zeros(size(Xuniq, 1), numel(genes));
    Gmat_count = zeros(1, numel(genes));
    
    matched_count = 0;
    matchedMask = false(1, numel(genes));
    
    if use_smart_mapping
        fprintf('    [Smart Mapping] Running with DIAGNOSTICS...\n');
        
        ignored_tokens = {'putative', 'fragment', 'complete', 'genome', 'reverse', ...
                          'forward', 'cds', 'rna', 'mrna', 'trna', 'rrna', 'orf', ...
                          'sequence', 'region', 'cluster', 'marker', 'variant', 'copy'};
        
        % Check if map is empty
        if geneIndexMap.Count == 0
            error('CRITICAL ERROR: geneIndexMap is empty! The model genes were not loaded.');
        end

        for k = 1:numel(uniqNames)
            dName = uniqNames(k);
            
            % Reverting to split on BOTH Hyphens and Underscores to be safe
            tokens = strsplit(dName, {'-', '_'}); 
            
            found_any = false;
            for t = 1:numel(tokens)
                tok = tokens{t};
                
                % Filter
                if strlength(tok) < 3 || ismember(lower(tok), ignored_tokens)
                    continue; 
                end

                key = lower(tok);
                if isKey(geneIndexMap, key)
                    idx = geneIndexMap(key);
                    
                    Gmat_acc(:, idx) = Gmat_acc(:, idx) + Xuniq(:, k);
                    Gmat_count(idx) = Gmat_count(idx) + 1;
                    
                    matchedMask(idx) = true;
                    found_any = true;
                    
                    
                    if matched_count == 0
                        fprintf('       [DEBUG SUCCESS] Mapped "%s" (from "%s") -> Model Gene Index %d\n', tok, dName, idx);
                    end
                else
                    
                    if k < 5
                         fprintf('       [DEBUG FAIL] Token "%s" (from "%s") not found in map.\n', tok, dName);
                    end
                end
            end
            if found_any, matched_count = matched_count + 1; end
        end
        
        % Averaging Step
        for g = 1:numel(genes)
            if Gmat_count(g) > 0
                Gmat(:, g) = Gmat_acc(:, g) ./ Gmat_count(g);
            end
        end
        
    else
        % Simple Mode
        for k = 1:numel(uniqNames)
            dName = uniqNames(k);
            if isKey(geneIndexMap, dName)
                idx = geneIndexMap(dName);
                Gmat(:, idx) = Xuniq(:, k);
                matchedMask(idx) = true;
                matched_count = matched_count + 1;
            end
        end
    end
    

    if sum(Gmat(:)) == 0
        warning('!!! OUTPUT MATRIX IS ALL ZEROS. NO GENES MATCHED !!!');
        fprintf('       Diagnostic: Checked %d input columns against %d model genes.\n', numel(uniqNames), geneIndexMap.Count);
    end

    % =====================================================================
    % SAVE
    % =====================================================================
    mappedDir = fullfile(dataDir, 'mapped_for_benchmarking');
    if ~exist(mappedDir, 'dir'), mkdir(mappedDir); end
    
    T_mapped = array2table(Gmat, 'VariableNames', cellstr(genes));
    mapped_outfile = fullfile(mappedDir, [baseName '_mapped_to_model.csv']);
    writetable(T_mapped, mapped_outfile);
    fprintf('    -> Saved Mapped Expression: %s\n', mapped_outfile);
    % =====================================================================

    Gmat = Gmat ./ gcount; 
    fprintf('Matched to model: %d input cols -> %d model genes (%.1f%% coverage)\n', ...
        matched_count, sum(matchedMask), sum(matchedMask)/numel(genes)*100);
    
    if matched_count == 0
        warning('Zero genes matched! Check use_smart_mapping flag or data headers.');
    end
    
    % E. Compute MRAS
    [mras, rxnNotBacked, rxnHasGPR] = compute_mras_batch_and_missing(postfixAll, Gmat, matchedMask);
    
    % F. Numerical Stability (Add microscopic noise to zero-variance)
    v = var(mras, 0, 1, 'omitnan');
    zero_var_mask = (v == 0) & rxnHasGPR;
    
    if any(zero_var_mask)
        noise_level = 1e-12 * max(std(mras, 0, 1));
        if noise_level == 0, noise_level = 1e-12; end
        for j = find(zero_var_mask)
            mras(:, j) = mras(:, j) + noise_level * abs(randn(size(mras, 1), 1));
        end
        fprintf('Added microscopic noise to %d zero-variance reactions\n', sum(zero_var_mask));
    end
    
    % G. Statistics & Output
    keep_for_cov = rxnHasGPR(:)';
    X_gpr_only   = mras(:, keep_for_cov);
    rxn_ids_gpr  = rxns(keep_for_cov);
    
    fprintf('Computing Covariance Matrix C...\n');
    cov_matrix = cov(X_gpr_only,1);
    
    fprintf('Computing Mean Outer Product (M_outer)...\n');
    mu_vec = mean(X_gpr_only, 1, 'omitnan')'; 
    mean_outer_matrix = mu_vec * mu_vec';    
    

    [n_cells, ~] = size(X_gpr_only);
    
    second_moment_matrix = (X_gpr_only' * X_gpr_only) / n_cells;

    if second_moment_mat
        fprintf('Check for Second moment Negativity: Smallest Value = %.4e | Count < 0 = %d\n', min(second_moment_matrix(:)), sum(second_moment_matrix(:) < 0));
    else
        secondmoment_temp=cov_matrix+mean_outer_matrix;
        fprintf('Check for Second moment Negativity: Smallest Value = %.4e | Count < 0 = %d\n', min(secondmoment_temp(:)), sum(secondmoment_temp(:) < 0));
        fprintf("There is a need for clamping before the optimisation!!!")
    end 

   
    nan_mras = sum(isnan(X_gpr_only), 'all');
    nan_cov  = sum(isnan(second_moment_matrix), 'all');
    
    if nan_mras > 0 || nan_cov > 0
        fprintf('!!! WARNING !!! NaN check failed for %s:\n', files(f).name);
        if nan_mras > 0, fprintf('    - MRAS Matrix (Reaction Activity): %d NaNs detected.\n', nan_mras); end
        if nan_cov > 0,  fprintf('    - Covariance Matrix: %d NaNs detected.\n', nan_cov); end
    else
        fprintf('NaN check passed: No missing values in reaction activity or covariance matrices.\n');
    end
    % ---------------------------------------------------------------------

    % --- SAVE OUTPUTS ---
    baseName = erase(files(f).name, '.csv');
   
    
    % 1. Mean Vector
    T_mean = table(mu_vec, 'VariableNames', {'MeanMRAS'});
    T_mean.Properties.RowNames = cellstr(rxn_ids_gpr);
    writetable(T_mean, fullfile(outDir, [baseName '_MEAN.csv']), 'WriteRowNames', true);
    
    second_moment = table(second_moment_matrix, 'VariableNames', {'SecondMoment'});
    second_moment.Properties.RowNames = cellstr(rxn_ids_gpr);
    writetable(second_moment, fullfile(outDir, [baseName '_SecondMoment.csv']), 'WriteRowNames', true);

    % 2. Covariance Matrix
    cov_table = array2table(cov_matrix, 'RowNames', cellstr(rxn_ids_gpr), 'VariableNames', cellstr(rxn_ids_gpr));
    writetable(cov_table, fullfile(outDir, [baseName '_COV.csv']), 'WriteRowNames', true);
    
    % 3. Mean Outer Product
    mo_table = array2table(mean_outer_matrix, 'RowNames', cellstr(rxn_ids_gpr), 'VariableNames', cellstr(rxn_ids_gpr));
    writetable(mo_table, fullfile(outDir, [baseName '_MEAN_OUTER.csv']), 'WriteRowNames', true);
    
    % 4. Reaction List
    fid = fopen(fullfile(outDir, [baseName '_gpr_reaction_ids.txt']), 'w');
    fprintf(fid, '%s\n', rxn_ids_gpr{:});
    fclose(fid);
    
    % 5. Missing Reaction Report
    missIdx = find(rxnNotBacked);
    if ~isempty(missIdx)
        fid = fopen(fullfile(outDir, [baseName '_missing_rxns.txt']), 'w');
        for k = 1:numel(missIdx)
            j = missIdx(k);
            reason = 'NO_DATA'; 
            if ~rxnHasGPR(j), reason = 'NO_GPR'; end
            rule = strtrim(grRules(j)); 
            if rule == "", rule = "(no GPR)"; end
            fprintf(fid, '%s\t%s\t%s\r\n', rxns(j), reason, rule);
        end
        fclose(fid);
    end
    
    fprintf('Completed %s in %.1f seconds\n', files(f).name, toc);
end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================
function gcount = count_gene_promiscuity_text(grRules, genes)
    L = lower(genes); 
    gcount = zeros(1, numel(genes));
    for j = 1:numel(grRules)
        r = lower(grRules{j}); 
        if isempty(r), continue; end
        r = regexprep(r, '\<and\>', ' '); 
        r = regexprep(r, '\<or\>', ' '); 
        r = regexprep(r, '[\(\)]', ' ');
        toks = strsplit(strtrim(r)); 
        if isempty(toks), continue; end
        [tf, loc] = ismember(toks, L); 
        loc = unique(loc(tf)); 
        gcount(loc) = gcount(loc) + 1;
    end
end

function postfix = gpr_to_postfix(rule, geneIndexMap)
    r = rule; 
    r = regexprep(r, '\[|\]', ''); 
    r = regexprep(r, '\<AND\>', 'and', 'ignorecase'); 
    r = regexprep(r, '\<OR\>', 'or', 'ignorecase'); 
    r = regexprep(lower(r), '\<and\>', '&'); 
    r = regexprep(r, '\<or\>', '|'); 
    r = regexprep(r, '([\(\)\&\|])', ' $1 ');
    
    tokens = strsplit(strtrim(r)); 
    ops = {}; 
    opTop = 0; 
    postfix = {}; 
    prec = containers.Map({'|', '&'}, [1 2]);
    
    for i = 1:numel(tokens)
        t = strtrim(tokens{i}); 
        if isempty(t), continue; end
        
        if strcmp(t, '(')
            opTop = opTop + 1; 
            ops{opTop} = t;
        elseif strcmp(t, ')')
            while opTop > 0 && ~strcmp(ops{opTop}, '(')
                postfix{end+1} = ops{opTop}; 
                opTop = opTop - 1; 
            end
            if opTop > 0 && strcmp(ops{opTop}, '('), opTop = opTop - 1; end
        elseif any(strcmp(t, {'&', '|'}))
            while opTop > 0 && ~strcmp(ops{opTop}, '(') && prec(ops{opTop}) >= prec(t)
                postfix{end+1} = ops{opTop}; 
                opTop = opTop - 1; 
            end
            opTop = opTop + 1; 
            ops{opTop} = t;
        else
            gn = strrep(strrep(t, '''', ''), '"', ''); 
            key = lower(gn);
            if isKey(geneIndexMap, key)
                postfix{end+1} = {'g', geneIndexMap(key)}; 
            else
                postfix{end+1} = {'g', 0}; 
            end
        end
    end
    while opTop > 0, postfix{end+1} = ops{opTop}; opTop = opTop - 1; end
end

function [mras, rxnNotBacked, rxnHasGPR] = compute_mras_batch_and_missing(postfixAll, Gmat, matchedMask)
    ns = size(Gmat, 1); 
    nr = numel(postfixAll); 
    mras = zeros(ns, nr); 
    rxnNotBacked = false(1, nr); 
    rxnHasGPR = false(1, nr);
    
    for j = 1:nr
        pf = postfixAll{j};
        if isempty(pf)
            mras(:, j) = 0; 
            rxnHasGPR(j) = false; 
            rxnNotBacked(j) = true;
        else
            rxnHasGPR(j) = true; 
            idx = extract_gene_indices(pf); 
            idx = idx(idx > 0);
            
            hasMatchedGene = false; 
            if ~isempty(idx), hasMatchedGene = any(matchedMask(idx)); end
            
            rxnNotBacked(j) = ~hasMatchedGene;
            mras(:, j) = eval_postfix_vec(pf, Gmat);
        end
    end
end

function idx = extract_gene_indices(postfix)
    idx = []; 
    for i = 1:numel(postfix)
        tok = postfix{i}; 
        if iscell(tok), idx = [idx tok{2}]; end
    end
    idx = unique(idx);
end

function v = eval_postfix_vec(postfix, Gmat)
    stack = {};
    for i = 1:numel(postfix)
        tok = postfix{i};
        if iscell(tok)
            idx = tok{2}; 
            if idx == 0
                val = zeros(size(Gmat, 1), 1); 
            else
                val = Gmat(:, idx); 
                val(isnan(val)) = 0; 
            end
            stack{end+1} = val;
        elseif ischar(tok)
            if strcmp(tok, '|')
                b = stack{end}; stack(end) = []; 
                a = stack{end}; stack(end) = []; 
                stack{end+1} = a + b; % OR = Sum
            elseif strcmp(tok, '&')
                b = stack{end}; stack(end) = []; 
                a = stack{end}; stack(end) = []; 
                stack{end+1} = min(a, b); % AND = Min
            end
        end
    end
    v = stack{end};
end