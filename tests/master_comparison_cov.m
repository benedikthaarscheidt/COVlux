% MASTER_COMPARISON_ULTIMATE_ALL_PLOTS_EXTENDED
% =========================================================================
% 1. Runs COVLUX, iMAT, FASTCORE.
% 2. WIRING: Covariance Analysis (Var, CorrM, Coh, Retained Var%).
% 3. TOPOLOGY: Lost Components + DUAL CONNECTIVITY (LCC% & AlgConn).
% 4. INTEGRITY: AA% (Protein Synthesis) AND MTA% (Core Machinery) under
%    both MAXIMAL and MINIMAL media – printed side by side.
% 5. CONSISTENCY: Falsepositives (Low Exp / High Flux).
% 6. REPORTING: Full CSV and Summary Table (by Method).
% 7. VISUALIZATION: All 10 Figures with both media overlaid in key plots.
% 8. CORRELATION: Overlaid Scatter Plots PDF with media distinction.
% =========================================================================

% --- ROBUST PROJECT ROOT FINDER ---
if exist(fullfile(pwd, 'config', 'config.json'), 'file')
    projectRoot = pwd;
    fprintf('Project Root Detected (Current Dir): %s\n', projectRoot);
else
    currentSearchPath = fileparts(mfilename('fullpath'));
    projectRoot = '';
    while length(currentSearchPath) > 1
        if exist(fullfile(currentSearchPath, 'config', 'config.json'), 'file')
            projectRoot = currentSearchPath;
            fprintf('Project Root Detected (Relative to Script): %s\n', projectRoot);
            break;
        end
        newPath = fileparts(currentSearchPath);
        if strcmp(newPath, currentSearchPath), break; end
        currentSearchPath = newPath;
    end
end
if isempty(projectRoot)
    error(['Could not locate "config/config.json".\n' ...
           'Please run this script from the project root folder.']);
end

% --- LOAD CONFIG ---
configFile = fullfile(projectRoot, 'config', 'config.json');
config = jsondecode(fileread(configFile));

%% 1. PARAMETERS & PATHS
run_name_cluster = config.params.input_clustering_folder;
use_conditions   = isfield(config.params, 'use_conditions') && config.params.use_conditions;
targetBiomassName = 'BIOMASS_Ec_iML1515_WT_75p37M';
modelPath = fullfile(projectRoot, config.paths.models_dir, config.model.model_file);

if config.params.use_big_basis
    covluxBase = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_bigbasis');
else
    covluxBase = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_smallbasis');
end

% AUTO-DETECT LATEST RUN
d = dir(fullfile(covluxBase, 'Run_*'));
d = d([d.isdir]);
if isempty(d)
    error('No COVLUX run folders found in: %s', covluxBase);
end
[~, idx] = sort([d.datenum], 'descend');
latestRun = d(idx(1)).name;
resultsDir = fullfile(covluxBase, latestRun);
fprintf('Targeting COVLUX Results: %s\n', resultsDir);

% Covariance Input Directory (MRAS Outputs)
clusteringBase = fullfile(projectRoot, config.paths.results_dir, 'Clustering');
if use_conditions
    covInputDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition','mapped_for_benchmarking');
else
    covInputDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files','mapped_for_benchmarking');
end
if ~exist(covInputDir, 'dir')
    fprintf('Warning: MRAS output dir not found at %s. Checking parent...\n', covInputDir);
    covInputDir = fileparts(covInputDir);
end
fprintf('MRAS/Covariance Source:   %s\n', covInputDir);
fprintf('Expression Source:        %s\n', expressionDir);

% --- PHPP SETTINGS ---
n_points = 20;
glc_range = linspace(0.1, 20, n_points);
o2_range = linspace(0.1, 20, n_points);

% --- Media definitions (minimal substrates from heatmap script) ---
media_conditions = {'maximal', 'minimal'};
minimal_substrates = {'ex_glc__d_e_b', 'ex_glc_e_b', 'ex_o2_e_b', 'ex_h2o_e_b', 'ex_h_e_b', ...
                      'ex_nh4_e_b', 'ex_pi_e_b', 'ex_so4_e_b', 'ex_k_e_b', 'ex_na1_e_b', ...
                      'ex_mg2_e_b', 'ex_ca2_e_b', 'ex_cl_e_b', 'ex_fe2_e_b', 'ex_fe3_e_b', ...
                      'ex_zn2_e_b', 'ex_mn2_e_b', 'ex_cu2_e_b', 'ex_cobalt2_e_b', ...
                      'ex_mobd_e_b', 'ex_ni2_e_b', 'ex_sel_e_b', 'ex_tungs_e_b', ...
                      'ex_thm_e_b', 'ex_cbl1_e_b', 'ex_niacin_e_b', ...
                      'ex_glyc_e_b', 'ex_pyr_e_b', 'ex_ade_e_b', 'ex_hxan_e_b', 'ex_ura_e_b', ...
                      'ex_glu__L_e_b', 'ex_asp__L_e_b', 'ex_ser__L_e_b','ex_leu__L_e_b'};

minimal_substrates = {'ex_glc__d_e_b', 'ex_nh4_e_b', 'ex_pi_e_b', 'ex_so4_e_b', ...
                      'ex_o2_e_b', 'ex_h2o_e_b', 'ex_h_e_b', 'ex_k_e_b', ...
                      'ex_na1_e_b', 'ex_mg2_e_b', 'ex_ca2_e_b', 'ex_cl_e_b', ...
                      'ex_fe2_e_b', 'ex_fe3_e_b', 'ex_zn2_e_b', 'ex_mn2_e_b', ...
                      'ex_cu2_e_b', 'ex_cobalt2_e_b', 'ex_mobd_e_b', ...
                      'ex_ni2_e_b', 'ex_sel_e_b', 'ex_tungs_e_b'};
%% 2. LOAD MODEL & PREPARE SUBSYSTEMS
fprintf('Loading Base Model...\n');
m_struct = load(modelPath);
if isfield(m_struct, 'pruned_ir')
    model = m_struct.pruned_ir;
else
    fn = fieldnames(m_struct);
    model = m_struct.(fn{1});
end
[m, n] = size(model.S);
n_rxns_orig = n;
n_mets = m;
model_base = model;

% Run FBA to verify growth on original model
model_test = changeObjective(model_base, targetBiomassName);
sol = optimizeCbModel(model_test, 'max');
if sol.f > 1e-6
    fprintf('   [SUCCESS] The pruned model grows! Growth Rate: %.4f\n', sol.f);
else
    fprintf('   [FAILURE] The pruned model DOES NOT grow (Rate: %.4e).\n', sol.f);
    warning('Do not use this model. Check if essential reactions were pruned in Step 5.');
end

% --- SUBSYSTEM CLEANUP ---
if isfield(model_base, 'subSystems')
    raw_subs = model_base.subSystems;
    if ischar(raw_subs), raw_subs = cellstr(raw_subs); end
    flat_subs = strings(n_rxns_orig, 1);
    for i = 1:n_rxns_orig
        if isempty(raw_subs{i})
            flat_subs(i) = "Unassigned";
        elseif iscell(raw_subs{i})
            flat_subs(i) = string(raw_subs{i}{1});
        else
            flat_subs(i) = string(raw_subs{i});
        end
    end
else
    flat_subs = repmat("Unassigned", n_rxns_orig, 1);
end
u_subs = unique(flat_subs);
total_counts_map = containers.Map(cellstr(u_subs), zeros(length(u_subs),1));
for i = 1:length(u_subs)
    total_counts_map(char(u_subs(i))) = sum(flat_subs == u_subs(i));
end

% Setup Media & Objective
bio_idx = find(strcmp(model_base.rxns, targetBiomassName));
if isempty(bio_idx), bio_idx = find(contains(model_base.rxns, 'BIOMASS'),1); end
model_base.c = zeros(n_rxns_orig,1); model_base.c(bio_idx) = 1;

% Setup PhPP Indices
glc_idx = find(contains(model_base.rxns, 'EX_glc__D_e') & contains(model_base.rxns, '_b'), 1);
o2_idx  = find(contains(model_base.rxns, 'EX_o2_e') & contains(model_base.rxns, '_b'), 1);
do_phpp = ~isempty(glc_idx) && ~isempty(o2_idx);

% Create maximal and minimal base models
model_max = model_base;   % all exchanges open (default)

model_min = model_base;
exc_idx = find(contains(model_min.rxns, 'EX_'));   % Identify exchange reactions
model_min.lb(exc_idx) = 0;
model_min.ub(exc_idx) = 0;                         % Close all exchanges

% Open only the minimal substrates
for i = 1:length(minimal_substrates)
    rxn_id = find(strcmpi(model_min.rxns, minimal_substrates{i})); % case-insensitive
    if ~isempty(rxn_id)
        model_min.lb(rxn_id) = 0;
        model_min.ub(rxn_id) = 1000;
    else
        warning('Minimal substrate %s not found in model.', minimal_substrates{i});
    end
end
% Set specific uptake limits for glucose and oxygen
glc_idx_min = find(strcmpi(model_min.rxns, 'ex_glc__d_e_b') | strcmpi(model_min.rxns, 'ex_glc_e_b'));
if ~isempty(glc_idx_min)
    model_min.ub(glc_idx_min) = 10;   % 10 mmol/gDW·h
end
o2_idx_min = find(strcmpi(model_min.rxns, 'ex_o2_e_b'));
if ~isempty(o2_idx_min)
    model_min.ub(o2_idx_min) = 20;    % 20 mmol/gDW·h
end

% Define tasks for AA% and MTA% (used repeatedly)
aa_tasks = {'ala__L_c','arg__L_c','asn__L_c','asp__L_c','cys__L_c','gln__L_c','glu__L_c','gly_c',...
            'his__L_c','ile__L_c','leu__L_c','lys__L_c','met__L_c','phe__L_c','pro__L_c',...
            'ser__L_c','thr__L_c','trp__L_c','tyr__L_c','val__L_c'};
mta_tasks = {'atp_c', 'ctp_c', 'gtp_c', 'utp_c', 'datp_c', 'dctp_c', 'dgtp_c', 'dttp_c', ...
             'nad_c', 'nadh_c', 'nadp_c', 'nadph_c', 'fad_c', 'coa_c', 'pheme_c', 'pydx5p_c', ...
             'pyr_c', 'accoa_c', 'succoa_c', 'akg_c', 'oaa_c'};

%% 3. MAIN LOOP
files = dir(fullfile(resultsDir, '*_unique_lost_rxns_full.csv'));
stats_out = struct([]);
opts_lp = optimoptions('linprog','Display','none');

% --- HEADER (now includes both maximal and minimal columns) ---
fprintf('\n%-28s | %-5s | %-5s | %-8s | %-8s | %-6s | %-8s | %-6s | %-5s | %-6s | %-5s | %-5s | %-5s | %-5s | %-5s | %-8s | %-8s | %-5s | %-5s | %-6s | %-30s\n', ...
    'Cluster (Method)', 'Keep', 'Lost', 'Bio_max', 'Bio_min', 'LCC%', 'AlgConn', 'FVA%', 'Block', 'FP', 'PhPP', 'AA_max', 'AA_min', 'MTA_max', 'MTA_min', 'DeadE', 'Var(X)', 'CorrM', 'Coh.', 'nCC(L)', 'Sz(L)%', 'Most_Pruned_Cat');
fprintf('%s\n', repmat('-', 1, 300));

for k = 1:length(files)
    filename = files(k).name;
    cluster_full_name = strrep(filename, '_unique_lost_rxns_full.csv', '');
    expr_name_parts = strsplit(cluster_full_name, '_MRAS');
    expr_short_name = expr_name_parts{1};
    cluster_full_name = strrep(cluster_full_name, '_logCPM', '');

    % --- A. LOAD COVARIANCE ---
    covFile = fullfile(covInputDir, [cluster_full_name '_COV.csv']);
    if ~exist(covFile, 'file')
        d = dir(fullfile(covInputDir, [expr_short_name '*_COV*.csv']));
        if ~isempty(d), covFile = fullfile(d(1).folder, d(1).name); end
    end
    if exist(covFile, 'file')
        T_cov = readtable(covFile, 'ReadRowNames', true, 'PreserveVariableNames', true);
        rxnNames_X = string(T_cov.Properties.RowNames);
        X_full = table2array(T_cov);
        [~, map_X_to_Model] = ismember(rxnNames_X, model_base.rxns);
        has_covariance = true;
    else
        has_covariance = false;
    end

    % --- B. PREPARE EXPRESSION ---
    exprFile = fullfile(expressionDir, [expr_short_name, '_mapped_to_model.csv']);
    RH_indices = []; RL_indices = [];
    current_RL_set = [];
    if exist(exprFile, 'file')
        T = readtable(exprFile, 'VariableNamingRule', 'preserve');
        file_genes = T.Properties.VariableNames;
        gene_vals  = mean(T{:, :}, 1, 'omitnan')';
        valid_genes_struct = struct();
        valid_genes_struct.gene = file_genes(:);
        valid_genes_struct.value = gene_vals(:);
        [RH_indices, RL_indices] = getExpressionSets(model_base, valid_genes_struct);
        current_RL_set = RL_indices;
    else
        warning('Mapped file not found: %s. Please re-run MRAS script.', exprFile);
    end
    if ~isempty(bio_idx)
        RH_indices = unique([RH_indices; bio_idx]);
    end

    % --- C. GET SURVIVORS (independent of medium) ---
    survivors = cell(1, 3);
    try
        raw_lines = readlines(fullfile(files(k).folder, filename));
        if length(raw_lines)>0 && (contains(raw_lines(1), 'lost') || contains(raw_lines(1), 'Var')), raw_lines(1)=[]; end
        covlux_names = strtrim(raw_lines(strlength(raw_lines)>0));
        [~, lost_covlux] = ismember(covlux_names, model_base.rxns);
        lost_idx = lost_covlux(lost_covlux > 0);
        survivors{1} = setdiff(1:n_rxns_orig, lost_idx)';
    catch, survivors{1} = []; end

    try survivors{2} = run_imat_milp_strict(model_base, RH_indices, RL_indices); catch, survivors{2} = []; end
    try kept_fastcore = fastcore(model_base, RH_indices, 1e-4); [is_kept, ~] = ismember(model_base.rxns, kept_fastcore.rxns); survivors{3} = find(is_kept); catch, survivors{3} = []; end

    % --- D. ANALYZE EACH METHOD ---
    names = {'COVLUX', 'iMAT', 'FASTCORE'};
    for m_idx = 1:3
        curr_survivors_idx = survivors{m_idx};
        % Ensure biomass reaction is kept
        if ~ismember(bio_idx, curr_survivors_idx)
            curr_survivors_idx = union(curr_survivors_idx, bio_idx);
        end
        lost_indices = setdiff(1:n_rxns_orig, curr_survivors_idx);
        curr_name = names{m_idx};
        num_kept = length(curr_survivors_idx);
        num_lost = n_rxns_orig - num_kept;

        % ========== MEDIA-INDEPENDENT METRICS (computed once per method) ==========
        % --- 1. TOPOLOGY OF LOST ---
        num_lost_components = 0; largest_lost_fraction = 0;
        if num_lost > 0
            S_lost = model_max.S(:, lost_indices);
            connected_mets = any(S_lost ~= 0, 2);
            S_lost = S_lost(connected_mets, :);
            if ~isempty(S_lost) && size(S_lost, 2) > 0
                n_mets_lost = size(S_lost, 1);
                n_rxns_lost = size(S_lost, 2);
                Adj_lost = [sparse(n_mets_lost, n_mets_lost), abs(S_lost) ~= 0;
                           (abs(S_lost) ~= 0)', sparse(n_rxns_lost, n_rxns_lost)];
                if nnz(Adj_lost) == 0
                    num_lost_components = n_rxns_lost;
                    largest_lost_fraction = 100 / n_rxns_lost;
                else
                    [bins, ~] = conncomp(graph(Adj_lost));
                    rxn_bins = bins(n_mets_lost+1 : end);
                    unique_bins = unique(rxn_bins);
                    num_lost_components = length(unique_bins);
                    if num_lost_components > 0
                       [~,~,bin_idx] = unique(rxn_bins);
                       bin_counts = accumarray(bin_idx, 1);
                       largest_lost_fraction = (max(bin_counts) / n_rxns_lost) * 100;
                    end
                end
            else
                num_lost_components = num_lost;
                largest_lost_fraction = 100 / num_lost;
            end
        end

        % --- 2. SUBSYSTEM CONTEXT ---
        most_pruned_sub = "None"; most_pruned_pct = 0;
        if ~isempty(lost_indices)
            lost_subs = flat_subs(lost_indices);
            u_lost_subs = unique(lost_subs);
            max_pct = -1;
            for s = 1:length(u_lost_subs)
                this_sub = u_lost_subs(s);
                if this_sub == "Unassigned" || this_sub == "" || contains(this_sub, "Exchange", 'IgnoreCase', true)
                    continue;
                end
                count_lost = sum(lost_subs == this_sub);
                count_total = total_counts_map(char(this_sub));
                if count_total > 5
                    pct_lost = count_lost / count_total;
                    if pct_lost > max_pct
                        max_pct = pct_lost;
                        most_pruned_sub = this_sub;
                    end
                end
            end
            if max_pct > -1, most_pruned_pct = max_pct; end
        end

        % --- 3. WIRING METRICS ---
        var_per_rxn = 0; coherence = 0; corr_mass = 0;
        if has_covariance && num_kept > 1
            keep_X_mask = ismember(map_X_to_Model, curr_survivors_idx);
            if sum(keep_X_mask) > 1
                X_sub = X_full(keep_X_mask, keep_X_mask);
                var_total = trace(X_sub);
                var_per_rxn = var_total / size(X_sub, 1);
                off_diag = X_sub - diag(diag(X_sub));
                corr_mass = sum(abs(off_diag(:)));
                coherence = corr_mass / var_total;
            end
        end

        % --- 4. DUAL CONNECTIVITY (LCC% and AlgConn) ---
        S_kept = model_max.S;
        S_kept(:, lost_indices) = 0;
        Adj = [sparse(n_mets, n_mets), abs(S_kept) ~= 0;
               (abs(S_kept) ~= 0)', sparse(n_rxns_orig, n_rxns_orig)];
        G_kept = graph(Adj);
        bins = conncomp(G_kept);
        bin_counts = groupcounts(bins');
        [max_comp_size, big_bin_idx] = max(bin_counts);
        lcc_pct = (max_comp_size / (n_mets + n_rxns_orig)) * 100;

        alg_conn = 0;
        if lcc_pct > 0.1
            try
                nodes_lcc = find(bins == big_bin_idx);
                if length(nodes_lcc) > 2
                    G_sub = subgraph(G_kept, nodes_lcc);
                    L = laplacian(G_sub);
                    evals = eigs(L + 1e-9*speye(size(L)), 2, 'smallestabs');
                    evals = sort(real(evals));
                    if length(evals) >= 2, alg_conn = evals(2); end
                end
            catch, alg_conn = 0; end
        end

        % --- 5. FVA & FALSE POSITIVES (using maximal model) ---
        model_test_fva = model_max;
        if ~isempty(lost_indices)
            model_test_fva.lb(lost_indices) = 0;
            model_test_fva.ub(lost_indices) = 0;
        end
        blocked_count = 0; falsepositives_count = 0;
        if num_kept > 0
            tol_fva = 1e-6;
            for i = 1:num_kept
                rxn_id = curr_survivors_idx(i);
                c_fva = zeros(n_rxns_orig, 1);
                c_fva(rxn_id) = -1;
                [~, f_max, flag_fva] = linprog(c_fva, [], [], model_test_fva.S, zeros(n_mets,1), model_test_fva.lb, model_test_fva.ub, opts_lp);
                v_max_val = abs((flag_fva==1) * f_max);
                if v_max_val < tol_fva
                    blocked_count = blocked_count + 1;
                elseif ismember(rxn_id, current_RL_set)
                    falsepositives_count = falsepositives_count + 1;
                end
            end
            functional_ratio = (num_kept - blocked_count) / num_kept;
        else
            functional_ratio = 0;
        end

        % --- 6. DEAD-END METABOLITES (using structure only) ---
        has_prod = any(S_kept > 0, 2);
        has_cons = any(S_kept < 0, 2);
        internal_mets = ~contains(model_max.mets, '_e') & ~contains(model_max.mets, '_b');
        dead_ends = internal_mets & ((has_prod & ~has_cons) | (~has_prod & has_cons));
        num_dead_ends = sum(dead_ends);

        % ========== LOOP OVER MEDIA FOR BIOMASS, AA%, MTA% ==========
        % We'll compute for both media and store results
        for med_idx = 1:2   % 1=maximal, 2=minimal
            if med_idx == 1
                model_medium = model_max;
            else
                model_medium = model_min;
            end

            % Apply lost indices to this medium
            model_test = model_medium;
            if ~isempty(lost_indices)
                model_test.lb(lost_indices) = 0;
                model_test.ub(lost_indices) = 0;
            end

            % --- Biomass ---
            [~, fval, flag] = linprog(-model_test.c, [], [], model_test.S, zeros(n_mets,1), model_test.lb, model_test.ub, opts_lp);
            if flag == 1
                biomass_val = -fval;
            else
                biomass_val = 0;
            end
            if med_idx == 1
                biomass_max = biomass_val;
            else
                biomass_min = biomass_val;
            end

            % --- AA% (binary) ---
            passed_aa = 0;
            for ti = 1:length(aa_tasks)
                metIdx = find(strcmp(model_test.mets, aa_tasks{ti}));
                if ~isempty(metIdx)
                    model_check = model_test;
                    model_check.S = [model_check.S, sparse(n_mets, 1)];
                    model_check.S(metIdx, end) = -1;
                    model_check.c = [zeros(n_rxns_orig,1); 1];
                    model_check.lb = [model_check.lb; 0];
                    model_check.ub = [model_check.ub; 1000];
                    [~, f_aa, fl_aa] = linprog(-model_check.c, [], [], model_check.S, zeros(n_mets,1), model_check.lb, model_check.ub, opts_lp);
                    if fl_aa == 1 && -f_aa > 1e-3
                        passed_aa = passed_aa + 1;
                    end
                end
            end
            aa_pct = (passed_aa / length(aa_tasks)) * 100;
            if med_idx == 1
                aa_max = aa_pct;
            else
                aa_min = aa_pct;
            end

            % --- MTA% (binary) ---
            passed_mta = 0;
            for ti = 1:length(mta_tasks)
                metIdx = find(strcmp(model_test.mets, mta_tasks{ti}));
                if ~isempty(metIdx)
                    model_check = model_test;
                    model_check.S = [model_check.S, sparse(n_mets, 1)];
                    model_check.S(metIdx, end) = -1;
                    model_check.c = [zeros(n_rxns_orig,1); 1];
                    model_check.lb = [model_check.lb; 0];
                    model_check.ub = [model_check.ub; 1000];
                    [~, f_mta, fl_mta] = linprog(-model_check.c, [], [], model_check.S, zeros(n_mets,1), model_check.lb, model_check.ub, opts_lp);
                    if fl_mta == 1 && -f_mta > 1e-3
                        passed_mta = passed_mta + 1;
                    end
                end
            end
            mta_pct = (passed_mta / length(mta_tasks)) * 100;
            if med_idx == 1
                mta_max = mta_pct;
            else
                mta_min = mta_pct;
            end

            % --- PhPP (only for maximal, and only if biomass_max > 1e-6) ---
            if med_idx == 1 && do_phpp && biomass_max > 1e-6
                phpp_vals = [];
                for i_g = 1:n_points
                    for j_o = 1:n_points
                        model_p = model_test;
                        model_p.ub(glc_idx) = glc_range(i_g);
                        model_p.ub(o2_idx) = o2_range(j_o);
                        [~, f_p, fl_p] = linprog(-model_p.c, [], [], model_p.S, zeros(n_mets,1), model_p.lb, model_p.ub, opts_lp);
                        if fl_p == 1, phpp_vals(end+1) = -f_p; end
                    end
                end
                if ~isempty(phpp_vals)
                    avg_phpp_max = mean(phpp_vals(phpp_vals > 0.01));
                else
                    avg_phpp_max = 0;
                end
            end
        end % media loop

        % --- PRINT ROW (single row per method) ---
        sub_str = sprintf('%s (%.3f)', most_pruned_sub, most_pruned_pct);
        if strlength(sub_str) > 30, sub_str = extractBefore(sub_str, 30); end
        row_name = sprintf('%s (%s)', cluster_full_name, curr_name);
        if strlength(row_name) > 28, row_name = extractBefore(row_name, 28); end

        fprintf('%-28s | %-5d | %-5d | %-8.4f | %-8.4f | %-6.1f | %-8.4f | %-6.1f | %-5d | %-6d | %-5.4f | %-5.0f | %-5.0f | %-5.0f | %-5.0f | %-5d | %-8.2e | %-8.2e | %-5.2f | %-5d | %-6.1f | %-30s\n', ...
            row_name, num_kept, num_lost, biomass_max, biomass_min, lcc_pct, alg_conn, functional_ratio*100, ...
            blocked_count, falsepositives_count, avg_phpp_max, aa_max, aa_min, mta_max, mta_min, num_dead_ends, var_per_rxn, corr_mass, coherence, ...
            num_lost_components, largest_lost_fraction, sub_str);

        % --- SAVE TO STATS_OUT ---
        s_idx = length(stats_out)+1;
        stats_out(s_idx).Cluster = cluster_full_name;
        stats_out(s_idx).Method = curr_name;
        stats_out(s_idx).KeptCount = num_kept;
        stats_out(s_idx).LostCount = num_lost;
        stats_out(s_idx).Biomass_max = biomass_max;
        stats_out(s_idx).Biomass_min = biomass_min;
        stats_out(s_idx).LCC_Pct = lcc_pct;
        stats_out(s_idx).Alg_Connectivity = alg_conn;
        stats_out(s_idx).FunctionalPct_FVA = functional_ratio * 100;
        stats_out(s_idx).BlockedSurvivors = blocked_count;
        stats_out(s_idx).Falsepositives = falsepositives_count;
        stats_out(s_idx).PhPP_Avg = avg_phpp_max;
        stats_out(s_idx).AA_max = aa_max;
        stats_out(s_idx).AA_min = aa_min;
        stats_out(s_idx).MTA_max = mta_max;
        stats_out(s_idx).MTA_min = mta_min;
        stats_out(s_idx).DeadEnd_Metabolites = num_dead_ends;
        stats_out(s_idx).AvgVariance = var_per_rxn;
        stats_out(s_idx).CorrelationMass = corr_mass;
        stats_out(s_idx).CoherenceRatio = coherence;
        stats_out(s_idx).Lost_NumComponents = num_lost_components;
        stats_out(s_idx).Lost_LargestComponentPct = largest_lost_fraction;
        stats_out(s_idx).MostPrunedSubsystem = most_pruned_sub;
        stats_out(s_idx).MostPrunedScore = most_pruned_pct;
    end
end

% Save full results
writetable(struct2table(stats_out), fullfile(resultsDir, 'FULL_METHOD_COMPARISON_WITH_BOTH_MEDIA.csv'));
fprintf('\nDone. Saved to FULL_METHOD_COMPARISON_WITH_BOTH_MEDIA.csv\n');

% Summary by Method (averages of max and min)
T_results = struct2table(stats_out);
T_results.Failed_max = T_results.Biomass_max < 1e-6;
T_results.Failed_min = T_results.Biomass_min < 1e-6;

% Group statistics including the failure flags
summary = grpstats(T_results, 'Method', {'mean', 'std', 'sum'}, ...
    'DataVars', {'Biomass_max', 'Biomass_min', 'AA_max', 'AA_min', 'MTA_max', 'MTA_min', ...
                 'Falsepositives', 'LCC_Pct', 'Alg_Connectivity', 'Failed_max', 'Failed_min'});

prof_report = table();
prof_report.Method = summary.Method;
prof_report.Success_Rate_max = 100 * (1 - (summary.sum_Failed_max ./ summary.GroupCount));
prof_report.Success_Rate_min = 100 * (1 - (summary.sum_Failed_min ./ summary.GroupCount));
prof_report.Avg_AA_max = summary.mean_AA_max;
prof_report.Avg_AA_min = summary.mean_AA_min;
prof_report.Avg_MTA_max = summary.mean_MTA_max;
prof_report.Avg_MTA_min = summary.mean_MTA_min;
prof_report.Avg_Biomass_max = summary.mean_Biomass_max;
prof_report.Avg_Biomass_min = summary.mean_Biomass_min;
prof_report.Avg_Connectivity_LCC = summary.mean_LCC_Pct;
prof_report.Avg_Robustness_AlgConn = summary.mean_Alg_Connectivity;
prof_report.Avg_Falsepositives = summary.mean_Falsepositives;
fprintf('\n\n=== SUMMARY REPORT (by Method) ===\n');
disp(prof_report);
writetable(prof_report, fullfile(resultsDir, 'SUMMARY_REPORT_BY_METHOD.csv'));

%% 8. VISUALIZATION (ALL 10 FIGURES) - SMART MEDIA FILTERING
fprintf('\n=== Generating Visualization Plots (Smart Media Filtering) ===\n');

% =========================================================================
% 1. BUILD THE LONG TABLE (T_long)
% =========================================================================
T_results = struct2table(stats_out);
T_long = table();
for med = {'max', 'min'}
    med_name = med{1};
    T_med = T_results;
    if strcmp(med_name, 'max')
        T_med.Biomass = T_med.Biomass_max;
        T_med.AA      = T_med.AA_max;
        T_med.MTA     = T_med.MTA_max;
    else
        T_med.Biomass = T_med.Biomass_min;
        T_med.AA      = T_med.AA_min;
        T_med.MTA     = T_med.MTA_min;
    end
    T_med.Media = repmat({med_name}, height(T_med), 1);
    T_long = [T_long; T_med];
end
T_long.Media = categorical(T_long.Media);

% =========================================================================
% 2. PRE-CALCULATIONS (Done on T_long to ensure all columns exist)
% =========================================================================
T_long.TotalCount  = T_long.KeptCount + T_long.LostCount;
T_long.PruningRate = T_long.LostCount ./ T_long.TotalCount;
T_long.FragRatio   = T_long.Lost_NumComponents ./ (T_long.LostCount + eps);

% --- Yields for F-Scores ---
Yield_AA  = T_long.AA / 100;
Yield_MTA = T_long.MTA / 100;
Yield_LCC = T_long.LCC_Pct / 100;
Yield_Bio = T_long.Biomass ./ (max(T_long.Biomass) + eps);
Yield_Mod = max(0, 1 - T_long.FragRatio);

% --- F-Scores ---
T_long.F_AA  = 2 * (Yield_AA .* T_long.PruningRate)  ./ (Yield_AA + T_long.PruningRate + eps);
T_long.F_MTA = 2 * (Yield_MTA .* T_long.PruningRate) ./ (Yield_MTA + T_long.PruningRate + eps);
T_long.F_LCC = 2 * (Yield_LCC .* T_long.PruningRate) ./ (Yield_LCC + T_long.PruningRate + eps);
T_long.F_Bio = 2 * (Yield_Bio .* T_long.PruningRate) ./ (Yield_Bio + T_long.PruningRate + eps);
T_long.F_Mod = 2 * (Yield_Mod .* T_long.PruningRate) ./ (Yield_Mod + T_long.PruningRate + eps);

% --- Normalized Structural F-scores ---
maxVar = max(T_long.AvgVariance) + eps;
T_long.F_Var = 2 * ((T_long.AvgVariance ./ maxVar) .* T_long.PruningRate) ./ ((T_long.AvgVariance ./ maxVar) + T_long.PruningRate + eps);

maxCC = max(T_long.Lost_NumComponents) + eps;
T_long.F_nCC = 2 * ((1 - (T_long.Lost_NumComponents ./ maxCC)) .* T_long.PruningRate) ./ ((1 - (T_long.Lost_NumComponents ./ maxCC)) + T_long.PruningRate + eps);

% --- Efficiencies (Kept) ---
T_long.Eff_Kept_LCC = T_long.LCC_Pct ./ T_long.KeptCount;
T_long.Eff_Kept_AA  = T_long.AA ./ T_long.KeptCount;
T_long.Eff_Kept_MTA = T_long.MTA ./ T_long.KeptCount;
T_long.Eff_Kept_Mod = T_long.FragRatio ./ T_long.KeptCount;

% --- Efficiencies (Deleted) ---
T_long.Eff_Del_LCC = T_long.LCC_Pct ./ (T_long.LostCount + eps);
T_long.Eff_Del_AA  = T_long.AA ./ (T_long.LostCount + eps);
T_long.Eff_Del_MTA = T_long.MTA ./ (T_long.LostCount + eps);
T_long.Eff_Del_Mod = T_long.FragRatio ./ (T_long.LostCount + eps);

% =========================================================================
% 3. CREATE T_UNIQUE
% =========================================================================
% Extract only the 'max' rows. Because we do this AFTER all calculations, 
% T_unique safely inherits every single column (including .Method).
T_unique = T_long(T_long.Media == 'max', :);

% =========================================================================
% 4. CONSTANTS FOR PLOTTING
% =========================================================================
% List of functional variables that require the media split
media_dep_vars = {'AA', 'MTA', 'Biomass', 'PhPP_Avg', ...
                  'F_AA', 'F_MTA', 'F_Bio', ...
                  'Eff_Kept_AA', 'Eff_Kept_MTA', ...
                  'Eff_Del_AA', 'Eff_Del_MTA'};

% Colors and Labels
colors_method = [
    0      0.4470 0.7410;   % COVLUX (Blue)
    0.9290 0.6940 0.1250;   % FASTCORE (Yellow)
    0.8500 0.3250 0.0980    % iMAT (Red)
];
methods_list = unique(T_long.Method);
media_list   = categories(T_long.Media);

plotsDir = fullfile(resultsDir, 'plots');
if ~exist(plotsDir, 'dir'), mkdir(plotsDir); end
% --- FIGURE 1: TOPOLOGY BOXPLOTS ---
f1 = figure('Position', [100, 100, 1200, 800], 'Name', 'Topology Analysis', 'Color', 'w', 'Visible', 'off');
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; boxchart(categorical(T_long.Method), T_long.Alg_Connectivity, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('Alg. Connectivity (Kept)'); title('Model Robustness'); grid on; legend;
nexttile; boxchart(categorical(T_long.Method), T_long.LCC_Pct, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('LCC % (Kept)'); title('Network Integrity'); grid on;
nexttile; boxchart(categorical(T_long.Method), T_long.FragRatio, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('Fragmentation Ratio'); title('Deletion Modularity (Lower=Better)'); grid on;
nexttile; boxchart(categorical(T_long.Method), T_long.Lost_LargestComponentPct, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('Lost Chunk Size %'); title('Chunk Size'); grid on;
title(t, 'Fig 1: Topological Analysis (Media shown for completeness)');
saveas(f1, fullfile(plotsDir, 'Fig1_Topology_Boxplots.png'));
close(f1);

% =========================================================================
% FIGURE 2: METRICS VS DELETION SIZE (SCATTER)
% =========================================================================
f2 = figure('Position', [150, 150, 1400, 900], 'Name', 'Metrics vs Deletion Size', 'Color', 'w', 'Visible', 'off');
t2 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
metrics_to_plot = {'Alg_Connectivity', 'LCC_Pct', 'Lost_NumComponents', 'AA', 'Falsepositives', 'Biomass'};
titles_plot = {'Kept Connectivity', 'Kept Integrity LCC', 'Deleted Components', 'Amino Acid Yield ', 'False Positives', 'Biomass Prod.'};

for i = 1:length(metrics_to_plot)
    nexttile; hold on;
    is_dep = ismember(metrics_to_plot{i}, media_dep_vars);
    
    for m = 1:length(methods_list)
        if is_dep
            % Functional: Split by media
            for med = 1:length(media_list)
                idx = strcmp(T_long.Method, methods_list{m}) & (T_long.Media == media_list{med});
                if med == 1, m_shape = 'o'; else, m_shape = 's'; end
                scatter(T_long.LostCount(idx), T_long.(metrics_to_plot{i})(idx), 'SizeData', 40, 'Marker', m_shape, 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('%s (%s)', methods_list{m}, media_list{med}));
            end
        else
            % Structural: No media split
            idx = strcmp(T_unique.Method, methods_list{m});
            scatter(T_unique.LostCount(idx), T_unique.(metrics_to_plot{i})(idx), 'SizeData', 40, 'Marker', 'o', 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6, 'DisplayName', methods_list{m});
        end
    end
    xlabel('Number of DELETED Reactions'); ylabel(metrics_to_plot{i}, 'Interpreter', 'none'); title(titles_plot{i}); grid on;
    if i == 1, legend('Location', 'best', 'FontSize', 8); end
end
title(t2, 'Fig 2: Impact of Deletion Magnitude (○ maximal, □ minimal for functional)');
saveas(f2, fullfile(plotsDir, 'Fig2_Metrics_vs_DeletionSize.png')); close(f2);

% =========================================================================
% FIGURE 3: MODULARITY MECHANICS (Strictly Structural)
% =========================================================================
f3 = figure('Position', [150, 150, 1400, 600], 'Name', 'Modularity vs Deletion Size', 'Color', 'w', 'Visible', 'off');
t3 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: Frag Ratio
nexttile; hold on;
for m = 1:length(methods_list)
    idx = strcmp(T_unique.Method, methods_list{m});
    scatter(T_unique.LostCount(idx), T_unique.FragRatio(idx), 'SizeData', 60, 'Marker', 'o', 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6, 'DisplayName', methods_list{m});
end
xlabel('Deleted Reactions'); ylabel('Fragmentation Ratio'); title('Deletion Chaos vs Size'); grid on; legend('Location','best');

% Panel 2: Largest Chunk
nexttile; hold on;
for m = 1:length(methods_list)
    idx = strcmp(T_unique.Method, methods_list{m});
    scatter(T_unique.LostCount(idx), T_unique.Lost_LargestComponentPct(idx), 'SizeData', 60, 'Marker', 'o', 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6, 'DisplayName', methods_list{m});
end
xlabel('Deleted Reactions'); ylabel('Largest Chunk (%)'); title('Deletion Module Size'); grid on; legend('Location','best');
title(t3, 'Fig 3: Modularity Mechanics');
saveas(f3, fullfile(plotsDir, 'Fig3_Modularity_Mechanics.png')); close(f3);

% =========================================================================
% FIGURE 4: COMPREHENSIVE METRICS
% =========================================================================
f4 = figure('Position', [100, 100, 1600, 1000], 'Name', 'All Metrics vs Deleted Count', 'Color', 'w', 'Visible', 'off');
t4 = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
all_metrics = {'Biomass', 'AA', 'MTA', 'PhPP_Avg', 'Alg_Connectivity', 'LCC_Pct', 'Falsepositives', 'DeadEnd_Metabolites', 'AvgVariance', 'CorrelationMass', 'FragRatio', 'MostPrunedScore'};

for i = 1:length(all_metrics)
    nexttile; hold on;
    is_dep = ismember(all_metrics{i}, media_dep_vars);
    
    for m = 1:length(methods_list)
        if is_dep
            for med = 1:length(media_list)
                idx = strcmp(T_long.Method, methods_list{m}) & (T_long.Media == media_list{med});
                if med == 1, m_shape = 'o'; else, m_shape = 's'; end
                scatter(T_long.LostCount(idx), T_long.(all_metrics{i})(idx), 'SizeData', 30, 'Marker', m_shape, 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.5, 'DisplayName', sprintf('%s (%s)', methods_list{m}, media_list{med}));
            end
        else
            idx = strcmp(T_unique.Method, methods_list{m});
            scatter(T_unique.LostCount(idx), T_unique.(all_metrics{i})(idx), 'SizeData', 30, 'Marker', 'o', 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.5, 'DisplayName', methods_list{m});
        end
    end
    xlabel('Deleted Reactions'); ylabel(all_metrics{i}, 'Interpreter', 'none'); title(all_metrics{i}, 'Interpreter', 'none'); grid on;
    if i == 1, legend('Location','best','FontSize',6); end
end
title(t4, 'Fig 4: Comprehensive Metrics vs Deletion');
saveas(f4, fullfile(plotsDir, 'Fig4_AllMetrics_vs_Deletion.png')); close(f4);

% --- FIGURE 5: KEPT EFFICIENCY (Metric / Kept) ---
f5 = figure('Position', [100, 100, 1200, 800], 'Name', 'Kept Efficiency', 'Color', 'w', 'Visible', 'off');
t5 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; boxchart(categorical(T_long.Method), T_long.Eff_Kept_LCC, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('LCC / Kept'); title('Structural Density'); grid on;
nexttile; boxchart(categorical(T_long.Method), T_long.Eff_Kept_AA, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('AA / Kept'); title('Biosynthetic Density'); grid on;
nexttile; boxchart(categorical(T_long.Method), T_long.Eff_Kept_MTA, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('MTA / Kept'); title('Core Density'); grid on;
nexttile; boxchart(categorical(T_long.Method), T_long.FragRatio, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('Fragmentation Ratio'); title('Disorder (Reference)'); grid on;
title(t5, 'Fig 5: KEPT Efficiency (Function per Retained Unit) – Media comparison');
legend('Location','best');
saveas(f5, fullfile(plotsDir, 'Fig5_Kept_Efficiency_Boxplots.png'));
close(f5);

% --- FIGURE 6: DELETION EFFICIENCY (Metric / Deleted) ---
f6 = figure('Position', [100, 100, 1200, 800], 'Name', 'Deletion Ratios', 'Color', 'w', 'Visible', 'off');
t6 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; boxchart(categorical(T_long.Method), T_long.Eff_Del_LCC, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('LCC / Deleted'); title('Structural Retention vs Pruning'); grid on;
nexttile; boxchart(categorical(T_long.Method), T_long.Eff_Del_AA, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('AA Yield / Deleted'); title('Biosynthetic Retention vs Pruning'); grid on;
nexttile; boxchart(categorical(T_long.Method), T_long.Eff_Del_MTA, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('MTA Yield / Deleted'); title('Core Retention vs Pruning'); grid on;
nexttile; boxchart(categorical(T_long.Method), T_long.Eff_Del_Mod, 'GroupByColor', T_long.Media); colororder(lines(2)); ylabel('Fragmentation / Deleted'); title('Disorder per Deletion'); grid on;
title(t6, 'Fig 6: DELETION Ratios (Function per Pruned Unit) – Media comparison');
legend('Location','best');
saveas(f6, fullfile(plotsDir, 'Fig6_Deletion_Ratios_Boxplots.png'));
close(f6);

% =========================================================================
% FIGURE 7: DELETION SCALING
% =========================================================================
f7 = figure('Position', [100, 100, 1200, 800], 'Name', 'Deletion Scaling', 'Color', 'w', 'Visible', 'off');
t7 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
metrics_del = {'Eff_Del_LCC', 'Eff_Del_AA', 'Eff_Del_MTA', 'Eff_Del_Mod'};
titles_del = {'LCC / Deleted', 'AA / Deleted', 'MTA / Deleted', 'Fragmentation / Deleted'};

for i = 1:4
    nexttile; hold on;
    is_dep = ismember(metrics_del{i}, media_dep_vars);
    
    for m = 1:length(methods_list)
        if is_dep
            for med = 1:length(media_list)
                idx = strcmp(T_long.Method, methods_list{m}) & (T_long.Media == media_list{med});
                if med == 1, m_shape = 'o'; else, m_shape = 's'; end
                scatter(T_long.LostCount(idx), T_long.(metrics_del{i})(idx), 'SizeData', 50, 'Marker', m_shape, 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('%s (%s)', methods_list{m}, media_list{med}));
            end
        else
            idx = strcmp(T_unique.Method, methods_list{m});
            scatter(T_unique.LostCount(idx), T_unique.(metrics_del{i})(idx), 'SizeData', 50, 'Marker', 'o', 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6, 'DisplayName', methods_list{m});
        end
    end
    xlabel('Deleted Reactions'); ylabel(titles_del{i}, 'Interpreter', 'none'); title(titles_del{i}); grid on;
    if i == 1, legend('Location','best'); end
end
title(t7, 'Fig 7: Scaling of Deletion Efficiency ');
saveas(f7, fullfile(plotsDir, 'Fig7_Deletion_Scaling.png')); close(f7);

% =========================================================================
% FIGURE 8: KEPT EFFICIENCY SCALING
% =========================================================================
f8 = figure('Position', [100, 100, 1200, 800], 'Name', 'Kept Scaling', 'Color', 'w', 'Visible', 'off');
t8 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
metrics_kept = {'Eff_Kept_LCC', 'Eff_Kept_AA', 'Eff_Kept_MTA', 'Eff_Kept_Mod'};
titles_kept = {'LCC / Kept', 'AA / Kept', 'MTA / Kept', 'FragRatio / Kept'};

for i = 1:4
    nexttile; hold on;
    is_dep = ismember(metrics_kept{i}, media_dep_vars);
    
    for m = 1:length(methods_list)
        if is_dep
            for med = 1:length(media_list)
                idx = strcmp(T_long.Method, methods_list{m}) & (T_long.Media == media_list{med});
                if med == 1, m_shape = 'o'; else, m_shape = 's'; end
                scatter(T_long.LostCount(idx), T_long.(metrics_kept{i})(idx), 'SizeData', 50, 'Marker', m_shape, 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('%s (%s)', methods_list{m}, media_list{med}));
            end
        else
            idx = strcmp(T_unique.Method, methods_list{m});
            scatter(T_unique.LostCount(idx), T_unique.(metrics_kept{i})(idx), 'SizeData', 50, 'Marker', 'o', 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6, 'DisplayName', methods_list{m});
        end
    end
    xlabel('Deleted Reactions'); ylabel(titles_kept{i}, 'Interpreter', 'none'); title(titles_kept{i}); grid on;
    if i == 1, legend('Location','best'); end
end
title(t8, 'Fig 8: Scaling of Kept Efficiency)');
saveas(f8, fullfile(plotsDir, 'Fig8_Kept_Efficiency_Scaling.png')); close(f8);
methods_order = {'COVLUX', 'FASTCORE', 'iMAT'};

% =========================================================================
% FIGURE 9: BALANCED F-SCORES (4 METRICS)
% =========================================================================
f9 = figure('Position', [100, 100, 1200, 800], 'Name', 'Balanced F-Scores', 'Color', 'w', 'Visible', 'off');
t9 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 
f_metrics = {'F_Mod', 'F_AA', 'F_MTA', 'F_Bio'};
f_titles  = {'Modularity Balance', 'AA Yield Balance', 'MTA Yield Balance', 'Biomass Balance'};

for i = 1:4
    ax = nexttile; hold on;
    is_dep = ismember(f_metrics{i}, media_dep_vars);
    
    for m = 1:3
        m_name = methods_order{m};
        m_color = colors_method(m, :);
        
        if is_dep
            % Maximal Media (Solid, left shifted)
            idx_max = strcmp(T_long.Method, m_name) & (T_long.Media == 'max');
            if any(idx_max)
                bc = boxchart(m*ones(sum(idx_max),1) - 0.15, T_long.(f_metrics{i})(idx_max));
                set(bc, 'BoxFaceColor', m_color, 'MarkerColor', m_color, 'BoxFaceAlpha', 0.85, 'LineWidth', 1.2, 'BoxWidth', 0.25);
            end
            
            % Minimal Media (Transparent, right shifted)
            idx_min = strcmp(T_long.Method, m_name) & (T_long.Media == 'min');
            if any(idx_min)
                bc = boxchart(m*ones(sum(idx_min),1) + 0.15, T_long.(f_metrics{i})(idx_min));
                set(bc, 'BoxFaceColor', m_color, 'MarkerColor', m_color, 'BoxFaceAlpha', 0.25, 'LineWidth', 1.2, 'BoxWidth', 0.25);
            end
        else
            % Structural Metric (Solid, perfectly centered)
            idx = strcmp(T_unique.Method, m_name);
            if any(idx)
                bc = boxchart(m*ones(sum(idx),1), T_unique.(f_metrics{i})(idx));
                set(bc, 'BoxFaceColor', m_color, 'MarkerColor', m_color, 'BoxFaceAlpha', 0.85, 'LineWidth', 1.2, 'BoxWidth', 0.4);
            end
        end
    end
    
    % Format the manual X-axis
    xticks(1:3);
    xticklabels(methods_order);
    xlim([0.5, 3.5]);
    
    grid on; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
    ylabel('F1 Score', 'FontWeight', 'bold'); 
    title(f_titles{i}, 'FontSize', 12); 
end

% Create Custom Legend Patches (Neutral Gray to explain Opacity)
hold on;
hMax = patch(NaN, NaN, [0.4 0.4 0.4], 'FaceAlpha', 0.85, 'EdgeColor', 'k', 'LineWidth', 1.2);
hMin = patch(NaN, NaN, [0.4 0.4 0.4], 'FaceAlpha', 0.25, 'EdgeColor', 'k', 'LineWidth', 1.2);
lg = legend([hMax, hMin], {'Maximal Media', 'Minimal Media'}, 'Orientation', 'horizontal', 'FontSize', 11);
lg.Layout.Tile = 'north';
title(t9, 'Fig 9: Balanced Scores (Trade-off: Metric vs. Pruning)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(f9, fullfile(plotsDir, 'Fig9_Balanced_Scores.png')); close(f9);


% =========================================================================
% FIGURE 10: RAW METRICS (5 METRICS)
% =========================================================================
f10 = figure('Position', [100, 100, 1400, 800], 'Name', 'Raw Metrics Boxplots', 'Color', 'w', 'Visible', 'off');
t10 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
raw_metrics  = {'AvgVariance', 'Lost_NumComponents', 'FragRatio', 'AA', 'MTA'};
raw_y_labels = {'Trace(X)/N', 'Count', 'Fragmentation Ratio', 'Yield (%)', 'Yield (%)'};
raw_titles   = {'Captured Variance', 'Removed CCs', 'Modularity', 'Amino Acid Yield', 'Metabolic Task Yield'};

for i = 1:5
    ax = nexttile; hold on;
    is_dep = ismember(raw_metrics{i}, media_dep_vars);
    
    for m = 1:3
        m_name = methods_order{m};
        m_color = colors_method(m, :);
        
        if is_dep
            % Maximal Media
            idx_max = strcmp(T_long.Method, m_name) & (T_long.Media == 'max');
            if any(idx_max)
                bc = boxchart(m*ones(sum(idx_max),1) - 0.15, T_long.(raw_metrics{i})(idx_max));
                set(bc, 'BoxFaceColor', m_color, 'MarkerColor', m_color, 'BoxFaceAlpha', 0.85, 'LineWidth', 1.2, 'BoxWidth', 0.25);
            end
            
            % Minimal Media
            idx_min = strcmp(T_long.Method, m_name) & (T_long.Media == 'min');
            if any(idx_min)
                bc = boxchart(m*ones(sum(idx_min),1) + 0.15, T_long.(raw_metrics{i})(idx_min));
                set(bc, 'BoxFaceColor', m_color, 'MarkerColor', m_color, 'BoxFaceAlpha', 0.25, 'LineWidth', 1.2, 'BoxWidth', 0.25);
            end
        else
            % Structural Metric
            idx = strcmp(T_unique.Method, m_name);
            if any(idx)
                bc = boxchart(m*ones(sum(idx),1), T_unique.(raw_metrics{i})(idx));
                set(bc, 'BoxFaceColor', m_color, 'MarkerColor', m_color, 'BoxFaceAlpha', 0.85, 'LineWidth', 1.2, 'BoxWidth', 0.4);
            end
        end
    end
    
    % Format the manual X-axis
    xticks(1:3);
    xticklabels(methods_order);
    xlim([0.5, 3.5]);
    
    grid on; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
    ylabel(raw_y_labels{i}, 'FontWeight', 'bold'); 
    title(raw_titles{i}, 'FontSize', 12); 
end

% Create Custom Legend Patches
hold on;
hMax = patch(NaN, NaN, [0.4 0.4 0.4], 'FaceAlpha', 0.85, 'EdgeColor', 'k', 'LineWidth', 1.2);
hMin = patch(NaN, NaN, [0.4 0.4 0.4], 'FaceAlpha', 0.25, 'EdgeColor', 'k', 'LineWidth', 1.2);
lg = legend([hMax, hMin], {'Maximal Media', 'Minimal Media'}, 'Orientation', 'horizontal', 'FontSize', 11);
lg.Layout.Tile = 'north';
title(t10, 'Fig 10: Raw Metrics (Wiring, Topology, Modularity, Function)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(f10, fullfile(plotsDir, 'Fig10_Raw_Metrics.png')); close(f10);
%% 9. CORRELATION ANALYSIS (PDF with strict data parsing)
fprintf('=== Generating Correlation PDF Report ===\n');

metrics_list = {'Biomass', 'LCC_Pct', 'Alg_Connectivity', 'AA', 'MTA', 'Falsepositives', 'FragRatio', 'Eff_Del_AA'};
metric_names_clean = {'Biomass (Fun)', 'LCC (%) (Str)', 'Algebraic Connectivity (Str)', 'AA Yield (%) (Fun)', 'MTA Yield (%) (Fun)', 'False Positives (Str)', 'Fragmentation (Str)', 'Del. Efficiency AA (Fun)'};

pdfFile = fullfile(plotsDir, 'Method_Deletion_Correlations.pdf');
if exist(pdfFile, 'file'), delete(pdfFile); end

fig = figure('Visible', 'off', 'Position', [0 0 1200 1600], 'Color', 'w');
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
title(t, 'Correlation: Metrics vs DELETED Reactions', 'FontSize', 16, 'FontWeight', 'bold');

for i = 1:length(metrics_list)
    ax = nexttile; hold(ax, 'on');
    y_name = metrics_list{i};
    is_dep = ismember(y_name, media_dep_vars);
    
    legend_handles = [];
    legend_labels = {};
    
    for m = 1:length(methods_list)
        if is_dep
            % Functional: Split by media
            for med = 1:length(media_list)
                % STRICT string matching to handle 'Energy'/'Amino_Acids' instead of 'max'/'min'
                idx = strcmp(string(T_long.Method), string(methods_list{m})) & strcmp(string(T_long.Media), string(media_list{med}));
                sub_data = T_long(idx, :);
                
                if height(sub_data) > 0
                    % Force double and drop NaNs
                    x_vals = double(sub_data.LostCount(:));
                    y_vals = double(sub_data.(y_name)(:));
                    mask = ~isnan(x_vals) & ~isnan(y_vals);
                    x_vals = x_vals(mask);
                    y_vals = y_vals(mask);
                    
                    if ~isempty(x_vals)
                        if med == 1, m_shape = 'o'; else, m_shape = 's'; end
                        disp_name = sprintf('%s %s', methods_list{m}, media_list{med});
                        
                        h_sc = scatter(x_vals, y_vals, 'SizeData', 60, 'Marker', m_shape, 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6);
                        
                        if ~ismember(disp_name, legend_labels)
                            legend_labels{end+1} = disp_name;
                            legend_handles(end+1) = h_sc(1); 
                        end
                        
                        if length(x_vals) > 2 && std(y_vals) > 1e-6
                            [~, ~] = corr(x_vals, y_vals, 'Type', 'Spearman', 'Rows', 'complete');
                            p = polyfit(x_vals, y_vals, 1);
                            x_fit = linspace(min(x_vals), max(x_vals), 100);
                            y_fit = polyval(p, x_fit);
                            plot(x_fit, y_fit, '-', 'Color', colors_method(m,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
                        end
                    end
                end
            end
        else
            % Structural: No media split
            idx = strcmp(string(T_unique.Method), string(methods_list{m}));
            sub_data = T_unique(idx, :);
            
            if height(sub_data) > 0
                x_vals = double(sub_data.LostCount(:));
                y_vals = double(sub_data.(y_name)(:));
                mask = ~isnan(x_vals) & ~isnan(y_vals);
                x_vals = x_vals(mask);
                y_vals = y_vals(mask);
                
                if ~isempty(x_vals)
                    disp_name = char(methods_list{m});
                    h_sc = scatter(x_vals, y_vals, 'SizeData', 60, 'Marker', 'o', 'MarkerEdgeColor', colors_method(m,:), 'MarkerFaceColor', colors_method(m,:), 'MarkerFaceAlpha', 0.6);
                    
                    if ~ismember(disp_name, legend_labels)
                        legend_labels{end+1} = disp_name;
                        legend_handles(end+1) = h_sc(1);
                    end
                    
                    if length(x_vals) > 2 && std(y_vals) > 1e-6
                        [~, ~] = corr(x_vals, y_vals, 'Type', 'Spearman', 'Rows', 'complete');
                        p = polyfit(x_vals, y_vals, 1);
                        x_fit = linspace(min(x_vals), max(x_vals), 100);
                        y_fit = polyval(p, x_fit);
                        plot(x_fit, y_fit, '-', 'Color', colors_method(m,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    end
                end
            end
        end
    end
    
    xlabel('Number of DELETED Reactions'); ylabel(metric_names_clean{i}); grid on;
    
    % Only plot legend if we actually drew points
    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 7, 'Box', 'on');
    end
    
    title(metric_names_clean{i}, 'FontSize', 12); hold(ax, 'off');
end

exportgraphics(fig, pdfFile, 'ContentType', 'vector');
close(fig);
fprintf('Correlation PDF saved to %s\n', pdfFile);
%% FUNCTION: iMAT MILP (unchanged)
function [kept_indices] = run_imat_milp_strict(model, RH, RL)
    [m, n] = size(model.S);
    f = zeros(4*n, 1); f(n+RH)=1; f(3*n+RL)=1;
    c = -f;
    A_eq = [model.S, sparse(m, 3*n)]; b_eq = zeros(m,1);
    lb=[model.lb; zeros(3*n,1)]; ub=[model.ub; ones(3*n,1)];
    rows=[]; cols=[]; vals=[]; rhs=[]; cnt=0; eps=0.01;
    k=RH; nk=length(k);
    if nk>0, r=(cnt+1:cnt+nk)'; rows=[rows;r;r]; cols=[cols;k;n+k]; vals=[vals;-ones(nk,1);(eps-model.lb(k))]; rhs=[rhs;-model.lb(k)]; cnt=cnt+nk; end
    k=RL; nk=length(k);
    if nk>0, r=(cnt+1:cnt+nk)'; rows=[rows;r;r]; cols=[cols;k;3*n+k]; vals=[vals;ones(nk,1);model.ub(k)]; rhs=[rhs;model.ub(k)]; cnt=cnt+nk; end
    A_ineq = sparse(rows, cols, vals, cnt, 4*n);
    milp = struct('A',[A_eq; A_ineq], 'b',[b_eq; rhs], 'c',c, 'lb',lb, 'ub',ub, ...
        'csense',[repmat('E',m,1); repmat('L',cnt,1)], 'vartype',[repmat('C',n,1); repmat('B',3*n,1)], 'osense',1);
    try sol = solveCobraMILP(milp, 'timeLimit', 60); catch, sol = solveCobraMILP(milp); end
    if sol.stat==1, kept_indices=find(abs(sol.full(1:n))>1e-9); else, kept_indices=[]; end
end

function [RH_indices, RL_indices] = getExpressionSets(model, geneExpression)
    [rxnExpression, ~] = mapExpressionToReactions(model, geneExpression);
    has_gpr = rxnExpression > -1;
    gpr_indices = find(has_gpr);
    if isempty(gpr_indices)
        RH_indices = []; RL_indices = []; return;
    end
    gpr_expr = rxnExpression(gpr_indices);
    high_thresh = prctile(gpr_expr, 75);
    low_thresh = prctile(gpr_expr, 25);
    RH_tmp = gpr_indices(gpr_expr >= high_thresh);
    RL_tmp = gpr_indices(gpr_expr <= low_thresh);
    RH_indices = RH_tmp(:);
    RL_indices = RL_tmp(:);
    fprintf("Length of RH indices: %d \n", length(RH_indices));
end