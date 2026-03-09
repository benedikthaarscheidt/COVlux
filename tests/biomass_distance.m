%% COVLUX vs iMAT vs FASTCORE: GAP & VARIANCE ANALYSIS
% 1. GAP: Min reactions needed to restore Biomass.
% 2. VARIANCE: How much data variance is captured by the "Fixed" (Gap-filled) model?
% =========================================================================

% --- 1. ROBUST PATH FINDING ---
if exist(fullfile(pwd, 'config', 'config.json'), 'file')
    projectRoot = pwd;
else
    currentSearchPath = fileparts(mfilename('fullpath'));
    projectRoot = '';
    while length(currentSearchPath) > 1 
        if exist(fullfile(currentSearchPath, 'config', 'config.json'), 'file')
            projectRoot = currentSearchPath; break;
        end
        newPath = fileparts(currentSearchPath);
        if strcmp(newPath, currentSearchPath), break; end
        currentSearchPath = newPath;
    end
end
if isempty(projectRoot), error('Could not locate config.json'); end

% --- LOAD CONFIG ---
configFile = fullfile(projectRoot, 'config', 'config.json');
config = jsondecode(fileread(configFile));
targetBiomassName = 'BIOMASS_Ec_iML1515_WT_75p37M';

% --- PARAMETERS ---
usebigbasis      = config.params.use_big_basis;
usesecondmoment  = config.params.second_moment_mat;
run_name_cluster = config.params.input_clustering_folder;
use_conditions   = isfield(config.params, 'use_conditions') && config.params.use_conditions;

%% 2. PATH DEFINITIONS
if usebigbasis
    covluxBase = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_bigbasis');
else
    covluxBase = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_smallbasis');
end

clusteringBase = fullfile(projectRoot, config.paths.results_dir, 'Clustering');
if use_conditions
    covDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition');
else
    covDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files');
end

d = dir(fullfile(covluxBase, 'Run_*')); 
d = d([d.isdir]);
if isempty(d), error('No Run folders found in %s', covluxBase); end

folderNames = {d.name};
runTimes = NaT(length(d), 1);
for i = 1:length(d)
    timeStr = extractAfter(folderNames{i}, 'Run_');
    try
        runTimes(i) = datetime(timeStr, 'InputFormat', 'yyyy-MM-dd_HH-mm-ss');
    catch
        runTimes(i) = datetime(d(i).datenum, 'ConvertFrom', 'datenum');
    end
end
[~, sortIdx] = sort(runTimes, 'descend');
resultsDir = fullfile(covluxBase, d(sortIdx(1)).name);

if use_conditions
    covInputDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition','mapped_for_benchmarking');
else
    covInputDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files','mapped_for_benchmarking');
end
if ~exist(covDir, 'dir'), covDir = fileparts(covDir); end 

modelPath = fullfile(projectRoot, config.paths.models_dir, config.model.model_file);
fprintf('Results Source:    %s\n', resultsDir);
fprintf('Covariance Source: %s\n', covDir);

%% 3. LOAD MODEL
fprintf('Loading Model...\n');
data = load(modelPath);
if isfield(data, 'pruned_ir'), model = data.pruned_ir; else, vars=fieldnames(data); model=data.(vars{1}); end
n_rxns = length(model.rxns);

bio_idx = find(strcmp(model.rxns, targetBiomassName));
if isempty(bio_idx), bio_idx = find(contains(model.rxns, 'BIOMASS'),1); end
if isempty(bio_idx), error('Biomass reaction not found.'); end

%% MINIMAL MEDIA DEFINITION
minimal_substrates = {'ex_glc__d_e_b', 'ex_glc_e_b', 'ex_o2_e_b', 'ex_h2o_e_b', 'ex_h_e_b', ...
                      'ex_nh4_e_b', 'ex_pi_e_b', 'ex_so4_e_b', 'ex_k_e_b', 'ex_na1_e_b', ...
                      'ex_mg2_e_b', 'ex_ca2_e_b', 'ex_cl_e_b', 'ex_fe2_e_b', 'ex_fe3_e_b', ...
                      'ex_zn2_e_b', 'ex_mn2_e_b', 'ex_cu2_e_b', 'ex_cobalt2_e_b', ...
                      'ex_mobd_e_b', 'ex_ni2_e_b', 'ex_sel_e_b', 'ex_tungs_e_b', ...
                      'ex_thm_e_b', 'ex_cbl1_e_b', 'ex_niacin_e_b', ...
                      'ex_glyc_e_b', 'ex_pyr_e_b', 'ex_ade_e_b', 'ex_hxan_e_b', 'ex_ura_e_b', ...
                      'ex_glu__L_e_b', 'ex_asp__L_e_b', 'ex_ser__L_e_b','ex_leu__L_e_b'};
model_min = model;
exc_idx = find(startsWith(lower(model.rxns), 'ex_') & endsWith(lower(model.rxns), '_b'));
model_min.lb(exc_idx) = 0;
model_min.ub(exc_idx) = 0;      
for i = 1:length(minimal_substrates)
    rxn_id = find(strcmpi(model_min.rxns, minimal_substrates{i})); 
    if ~isempty(rxn_id)
        model_min.lb(rxn_id) = 0;
        model_min.ub(rxn_id) = 1000;
    else
        warning('Minimal substrate %s not found in model.', minimal_substrates{i});
    end
end

fprintf('Printing non-zero S-matrix coefficients for %d exchange reactions:\n', length(exc_idx));
fprintf('%-25s | %-20s | %-10s\n', 'Reaction ID', 'Metabolite', 'Coefficient');
fprintf('%s\n', repmat('-', 1, 60));

fprintf('\n--- Running Pre-Flight Check on Unpruned Model ---\n');
[~, dead_in_base] = perform_biomass_autopsy(model, 1:n_rxns, bio_idx);
bio_col = model.S(:, bio_idx);
all_precursors = model.mets(bio_col < -1e-6);
valid_precursors = setdiff(all_precursors, dead_in_base);
fprintf('   Base model can synthesize %d/%d precursors.\n', length(valid_precursors), length(all_precursors));
fprintf('   Skipping %d impossible cofactors/precursors in gap analysis.\n', length(dead_in_base));

for i = 1:length(exc_idx)
    rxn_idx = exc_idx(i);
    rxn_name = model_min.rxns{rxn_idx};
    met_indices = find(model_min.S(:, rxn_idx));
    for j = 1:length(met_indices)
        met_idx = met_indices(j);
        met_name = model_min.mets{met_idx};
        coeff = model_min.S(met_idx, rxn_idx);
        fprintf('%-25s | %-20s | %-10.2f\n', rxn_name, met_name, coeff);
    end
end 

%% 4. MAIN LOOP
files = dir(fullfile(resultsDir, '*_unique_lost_rxns_full.csv'));
if isempty(files), error('No result CSVs found.'); end
Stats = struct([]);
fprintf('\n%-20s | %-5s %-5s %-5s | %-5s %-5s %-5s | %-5s %-5s %-5s\n', 'Cluster', 'G_COVlux', 'G_FastCore', 'G_iMAT', 'V_COVlux', 'V_FastCore', 'V_iMAT', 'L_COVlux', 'L_FastCore', 'L_iMAT');
fprintf('%s\n', repmat('-',1,105));

for k = 1:length(files)
    filename = files(k).name;
    cluster_full_name = strrep(filename, '_unique_lost_rxns_full.csv', '');
    parts = strsplit(cluster_full_name, '_MRAS');
    short_name = parts{1};
    expr_name_parts = strsplit(cluster_full_name, '_MRAS');
    expr_short_name = expr_name_parts{1};
    
    % --- A. GET COVLUX SURVIVORS ---
    try
        raw_lines = readlines(fullfile(files(k).folder, filename));
        if length(raw_lines) > 0 && (contains(raw_lines(1), 'lost') || contains(raw_lines(1), 'Var'))
            raw_lines(1) = [];
        end
        covlux_names = strtrim(raw_lines(strlength(raw_lines) > 0));
        [~, lost_covlux] = ismember(covlux_names, model.rxns);
        lost_idx = lost_covlux(lost_covlux > 0);
        survivors_cov = setdiff(1:n_rxns, lost_idx)';
        lost_cov_count = length(lost_idx);
        kept_cov = false(n_rxns, 1);
        kept_cov(survivors_cov) = true;
    catch ME
        fprintf('Error reading COVlux file for %s: %s\n', short_name, ME.message);
        kept_cov = true(n_rxns, 1); 
        lost_cov_count = 0;
    end
    
    % --- B. GET iMAT/FASTCORE SURVIVORS ---
    exprFile = fullfile(expressionDir,  [expr_short_name, '_logCPM_mapped_to_model.csv']);
    RH = []; RL = [];
    if exist(exprFile, 'file')
        T = readtable(exprFile, 'VariableNamingRule', 'preserve');
        gene_vals = mean(T{:, :}, 1, 'omitnan')';
        valid_genes_struct = struct();
        valid_genes_struct.gene = T.Properties.VariableNames(:); 
        valid_genes_struct.value = gene_vals(:);
        [RH, RL] = getExpressionSets(model, valid_genes_struct);
    else
        fprintf('Warning: Mapped expression file not found for %s\n', short_name);
    end
    
    idx_imat = run_imat_milp_standard(model, RH, RL);
    kept_imat = false(n_rxns, 1); kept_imat(idx_imat) = true;
    lost_imat_count = n_rxns - sum(kept_imat);
    
    try
        fast_model = fastcore(model, RH, 1e-4);
        [is_kept, ~] = ismember(model.rxns, fast_model.rxns);
        kept_fast = is_kept;
        lost_fast_count = n_rxns - sum(kept_fast);
    catch
        kept_fast = false(n_rxns, 1);
        lost_fast_count = n_rxns;
    end
    
    % --- C. GAP ANALYSIS & AUTOPSY ---
    
    % 1. Calculate Functional Gaps (Maximal Media)
    [gap_cov, added_cov]   = calculate_gap_and_patch(model, kept_cov, bio_idx);
    [gap_imat, added_imat] = calculate_gap_and_patch(model, kept_imat, bio_idx);
    [gap_fast, added_fast] = calculate_gap_and_patch(model, kept_fast, bio_idx);
    methods_added = {added_cov, added_imat, added_fast};
    method_names  = {'COVLUX', 'iMAT', 'FASTCORE'};
    
    fprintf('\n   --- ADDED REACTIONS (Gap Fixes - Max Media) ---\n');
    for m_idx = 1:3
        curr_added = methods_added{m_idx};
        fprintf('   [%s] Gap: %d\n', method_names{m_idx}, length(curr_added));
        if ~isempty(curr_added)
            for r = 1:min(10, length(curr_added))
                rxn_id = curr_added(r);
                rxn_name = model.rxns{rxn_id};
                if isfield(model, 'subSystems')
                    sub = model.subSystems{rxn_id};
                    if iscell(sub), sub = sub{1}; end
                else
                    sub = 'N/A';
                end
                fprintf('      + %-15s | %s\n', rxn_name, sub);
            end
            if length(curr_added) > 10, fprintf('      ... and %d more.\n', length(curr_added) - 10); end
        else
            fprintf('      (No reactions added - network is self-sufficient)\n');
        end
    end
    fprintf('   -----------------------------------\n');
    
    % 2. Calculate Functional Gaps (Minimal Media)
    [gap_cov_min, added_cov_min]   = calculate_gap_and_patch(model_min, kept_cov, bio_idx);
    [gap_imat_min, added_imat_min] = calculate_gap_and_patch(model_min, kept_imat, bio_idx);
    [gap_fast_min, added_fast_min] = calculate_gap_and_patch(model_min, kept_fast, bio_idx);
    
    % 3. AUTOPSY
    idx_cov  = find(kept_cov);
    idx_imat = find(kept_imat);
    idx_fast = find(kept_fast);
    
    fprintf('\n   --- BIOMASS AUTOPSY (Max Media) ---\n');
    [intact_cov, dead_cov] = perform_biomass_autopsy(model, idx_cov, bio_idx);
    tot_pre = length(intact_cov) + length(dead_cov);
    fprintf('   [COVLUX] Intact: %2d/%2d | Dead: %2d\n', length(intact_cov), tot_pre, length(dead_cov));
    if ~isempty(dead_cov), fprintf('      Missing: %s\n', strjoin(dead_cov(1:min(5, end)), ', ')); end
    
    [intact_imat, dead_imat] = perform_biomass_autopsy(model, idx_imat, bio_idx);
    fprintf('   [iMAT]   Intact: %2d/%2d | Dead: %2d\n', length(intact_imat), tot_pre, length(dead_imat));
    if ~isempty(dead_imat), fprintf('      Missing: %s\n', strjoin(dead_imat(1:min(5, end)), ', ')); end
    
    [intact_fast, dead_fast] = perform_biomass_autopsy(model, idx_fast, bio_idx);
    fprintf('   [FAST]   Intact: %2d/%2d | Dead: %2d\n', length(intact_fast), tot_pre, length(dead_fast));
    if ~isempty(dead_fast), fprintf('      Missing: %s\n', strjoin(dead_fast(1:min(5, end)), ', ')); end
    fprintf('   ---------------------------------------\n');
    
    fprintf('\n   --- BIOMASS AUTOPSY (Minimal Media) ---\n');
    [intact_cov_min, dead_cov_min] = perform_biomass_autopsy(model_min, idx_cov, bio_idx);
    fprintf('   [COVLUX] Intact: %2d/%2d | Dead: %2d\n', length(intact_cov_min), tot_pre, length(dead_cov_min));
    if ~isempty(dead_cov_min), fprintf('      Missing: %s\n', strjoin(dead_cov_min(1:min(5, end)), ', ')); end
    
    [intact_imat_min, dead_imat_min] = perform_biomass_autopsy(model_min, idx_imat, bio_idx);
    fprintf('   [iMAT]   Intact: %2d/%2d | Dead: %2d\n', length(intact_imat_min), tot_pre, length(dead_imat_min));
    if ~isempty(dead_imat_min), fprintf('      Missing: %s\n', strjoin(dead_imat_min(1:min(5, end)), ', ')); end
    
    [intact_fast_min, dead_fast_min] = perform_biomass_autopsy(model_min, idx_fast, bio_idx);
    fprintf('   [FAST]   Intact: %2d/%2d | Dead: %2d\n', length(intact_fast_min), tot_pre, length(dead_fast_min));
    if ~isempty(dead_fast_min), fprintf('      Missing: %s\n', strjoin(dead_fast_min(1:min(5, end)), ', ')); end
    fprintf('   ---------------------------------------\n');
    
    % --- DUAL RECONSTITUTION VERIFICATION (FBA) ---
    fprintf('\n   --- RECONSTITUTION VERIFICATION (FBA) ---\n');
    
    % Infrastructure Mask for Safety
    rxns_lower = lower(model.rxns);
    infrastructure_mask = contains(rxns_lower, 'biomass') | ...
                          startsWith(rxns_lower, 'ex_') | ...
                          contains(rxns_lower, 'atpm');
    infra_idx = find(infrastructure_mask);
    methods_kept_idx = {idx_cov, idx_imat, idx_fast};

    % --- SCENARIO 1: MAXIMAL MEDIA TEST ---
    fprintf('   >> SCENARIO 1: MAXIMAL MEDIA (Rich Environment)\n');
    methods_added_max = {added_cov, added_imat, added_fast};
    for m_idx = 1:3
        check_model_max = model; 
        
        final_kept_max = unique([methods_kept_idx{m_idx}; methods_added_max{m_idx}; bio_idx]);
        
        is_lost = true(n_rxns, 1);
        is_lost(final_kept_max) = false;
        check_model_max.lb(is_lost) = 0;
        check_model_max.ub(is_lost) = 0;
        
        % FIX: Explicitly tell FBA to maximize Biomass and NOTHING ELSE
        check_model_max.c = zeros(n_rxns, 1);
        check_model_max.c(bio_idx) = 1;
        
        sol_max = optimizeCbModel(check_model_max, 'max');
        if sol_max.stat == 1 && sol_max.f > 1e-6
            fprintf('      [%s] VERIFIED: Grows at %.4f h^-1\n', method_names{m_idx}, sol_max.f);
        else
            fprintf('      [%s] FAILED: Flux %.2e (Solver Status: %d)\n', method_names{m_idx}, sol_max.f, sol_max.stat);
            
            % --- THE DIAGNOSTIC ---
            fprintf('         > Running FBA Diagnostic...\n');
            try
                relaxOpt = struct('internalRelax', 2, 'exchangeRelax', 2, 'bRelax', 0);
                [~, ~, ~, ~, ~, relaxDir] = relaxFBA(check_model_max, relaxOpt);
                bad_rxns = check_model_max.rxns(relaxDir > 1e-6);
                if ~isempty(bad_rxns)
                    fprintf('         > FBA starved for: %s\n', strjoin(bad_rxns(1:min(5, end)), ', '));
                else
                    fprintf('         > FBA failed due to an unbounded thermodynamic loop.\n');
                end
            catch
                fprintf('         > Diagnostic failed to isolate the bottleneck.\n');
            end
        end
    end

    % --- SCENARIO 2: MINIMAL MEDIA TEST ---
    fprintf('\n   >> SCENARIO 2: MINIMAL MEDIA (Strict Substrates)\n');
    methods_added_min = {added_cov_min, added_imat_min, added_fast_min};
    for m_idx = 1:3
        check_model_min = model_min; 
        
        final_kept_min = unique([methods_kept_idx{m_idx}; methods_added_min{m_idx}; bio_idx]);
        
        is_lost = true(n_rxns, 1);
        is_lost(final_kept_min) = false;
        check_model_min.lb(is_lost) = 0;
        check_model_min.ub(is_lost) = 0;
        
        % FIX: Explicitly tell FBA to maximize Biomass and NOTHING ELSE
        check_model_min.c = zeros(n_rxns, 1);
        check_model_min.c(bio_idx) = 1;
        
        sol_min = optimizeCbModel(check_model_min, 'max');
        if sol_min.stat == 1 && sol_min.f > 1e-6
            fprintf('      [%s] VERIFIED: Grows at %.4f h^-1\n', method_names{m_idx}, sol_min.f);
        else
            fprintf('      [%s] FAILED: Flux %.2e (Solver Status: %d)\n', method_names{m_idx}, sol_min.f, sol_min.stat);
            
            % --- THE DIAGNOSTIC ---
            fprintf('         > Running FBA Diagnostic...\n');
            try
                relaxOpt = struct('internalRelax', 2, 'exchangeRelax', 2, 'bRelax', 0);
                [~, ~, ~, ~, ~, relaxDir] = relaxFBA(check_model_min, relaxOpt);
                bad_rxns = check_model_min.rxns(relaxDir > 1e-6);
                if ~isempty(bad_rxns)
                    fprintf('         > FBA starved for: %s\n', strjoin(bad_rxns(1:min(5, end)), ', '));
                else
                    fprintf('         > FBA failed due to an unbounded thermodynamic loop.\n');
                end
            catch
                fprintf('         > Diagnostic failed to isolate the bottleneck.\n');
            end
        end
    end
    fprintf('   -------------------------------------------\n');
    
    % --- D. VARIANCE ANALYSIS (On Fixed Models - Max Media) ---
    if usesecondmoment
        covFile = fullfile(covDir, [cluster_full_name '_logCPM_SecondMoment.csv']);
    else 
        covFile = fullfile(covDir, [cluster_full_name '_logCPM_COV.csv']);
    end 
    if ~exist(covFile, 'file'), covFile = fullfile(covDir, [short_name '_MRAS_COV.csv']); end
    
    var_cov = NaN; var_imat = NaN; var_fast = NaN;
    avg_cov = NaN; avg_imat = NaN; avg_fast = NaN;
    coh_cov = NaN; coh_imat = NaN; coh_fast = NaN; 
    conn_cov = NaN; conn_imat = NaN; conn_fast = NaN; 
    
    if exist(covFile, 'file')
        try
            T_cv = readtable(covFile, 'ReadRowNames', true, 'PreserveVariableNames', true);
            X = table2array(T_cv);
            rxns_X = string(T_cv.Properties.RowNames);
            [~, map_Model_to_X] = ismember(model.rxns, rxns_X);
            total_var = trace(X);
            std_vec = sqrt(diag(X));
            Corr_Matrix = X ./ (std_vec * std_vec');
            Corr_Matrix(isnan(Corr_Matrix)) = 0; 
            
            get_raw_sum = @(mask) sum(diag(X(map_Model_to_X(mask & map_Model_to_X>0), map_Model_to_X(mask & map_Model_to_X>0))));
            calc_pct = @(raw_sum) (raw_sum / total_var) * 100;
            get_coherence = @(mask) mean(abs(Corr_Matrix(map_Model_to_X(mask & map_Model_to_X>0), map_Model_to_X(mask & map_Model_to_X>0))), 'all', 'omitnan');
            get_connectivity = @(mask) mean(any(abs(Corr_Matrix(map_Model_to_X(mask & map_Model_to_X>0), map_Model_to_X(mask & map_Model_to_X>0))) > 0.5, 2));
            
            fixed_cov = kept_cov; fixed_cov(added_cov) = true; fixed_cov(added_cov_min) = true;
            raw_cov = get_raw_sum(fixed_cov); var_cov = calc_pct(raw_cov); 
            avg_cov = raw_cov / sum(fixed_cov); coh_cov = get_coherence(fixed_cov); conn_cov = get_connectivity(fixed_cov); 
            
            fixed_imat = kept_imat; fixed_imat(added_imat) = true;fixed_imat(added_imat_min) = true;
            raw_imat = get_raw_sum(fixed_imat); var_imat = calc_pct(raw_imat);
            avg_imat = raw_imat / sum(fixed_imat); coh_imat = get_coherence(fixed_imat); conn_imat = get_connectivity(fixed_imat);
            
            fixed_fast = kept_fast; fixed_fast(added_fast) = true;fixed_fast(added_fast_min) = true;
            raw_fast = get_raw_sum(fixed_fast); var_fast = calc_pct(raw_fast);
            avg_fast = raw_fast / sum(fixed_fast); coh_fast = get_coherence(fixed_fast); conn_fast = get_connectivity(fixed_fast); 
        catch ME
            warning('Metric calculation failed for %s: %s', short_name, ME.message);
        end
    end
    
    fprintf('%-20s | %-5d %-5d %-5d | %-5d %-5d %-5d  | %-5.1f %-5.1f %-5.1f | %-5d %-5d %-5d\n', ...
        short_name, gap_cov, gap_fast, gap_imat,gap_cov_min, gap_fast_min, gap_imat_min, var_cov, var_fast, var_imat, ...
        lost_cov_count, lost_fast_count, lost_imat_count);
    if ~isnan(coh_cov)
        fprintf('   > Coherence: COV=%.3f | iMAT=%.3f | FAST=%.3f\n', coh_cov, coh_imat, coh_fast);
    end
    
    Stats(k).Cluster = string(short_name);
    Stats(k).Gap_COV = gap_cov; Stats(k).Gap_FAST = gap_fast; Stats(k).Gap_iMAT = gap_imat;
    Stats(k).Gap_Min_COV = gap_cov_min; Stats(k).Gap_Min_iMAT = gap_imat_min; Stats(k).Gap_Min_FAST = gap_fast_min;
    Stats(k).Var_COV = var_cov; Stats(k).Var_FAST = var_fast; Stats(k).Var_iMAT = var_imat;
    Stats(k).AvgVar_COV = avg_cov; Stats(k).AvgVar_FAST = avg_fast; Stats(k).AvgVar_iMAT = avg_imat; 
    Stats(k).Lost_COV = lost_cov_count; Stats(k).Lost_FAST = lost_fast_count; Stats(k).Lost_iMAT = lost_imat_count;
    Stats(k).Coherence_COV = coh_cov; Stats(k).Coherence_FAST = coh_fast; Stats(k).Coherence_iMAT = coh_imat;
    Stats(k).Conn_COV = conn_cov; Stats(k).Conn_FAST = conn_fast; Stats(k).Conn_iMAT = conn_imat;
end

%% 5. VISUALIZATION
cluster_full_name = strrep(cluster_full_name, '_logCPM', '');
plotdir = fullfile(resultsDir, "plots");
if ~exist(plotdir, 'dir'), mkdir(plotdir); end
outputPdfPath = fullfile(plotdir, 'Functional_Gap_and_Variance.pdf');
T = struct2table(Stats);

fig = figure('Name', 'Gap vs Variance', 'Color', 'w', 'Position', [100 100 1200 1800]);
t = tiledlayout(5, 1, 'Padding', 'compact');

nexttile;
bar_data_gap = [T.Gap_COV, T.Gap_FAST, T.Gap_iMAT];
b1 = bar(bar_data_gap, 'grouped');
b1(1).FaceColor = [0 0.45 0.74]; b1(2).FaceColor = [0.93 0.69 0.13]; b1(3).FaceColor = [0.85 0.33 0.1];
ylabel('Reactions Added'); title('1. Max Media Gap (Distance to Biomass)');
legend({'COVlux', 'FASTCORE', 'iMAT'}, 'Location', 'best'); xticklabels(T.Cluster); xtickangle(45); grid on;

nexttile;
bar_data_min = [T.Gap_Min_COV, T.Gap_Min_FAST, T.Gap_Min_iMAT];
b5 = bar(bar_data_min, 'grouped');
b5(1).FaceColor = [0 0.45 0.74]; b5(2).FaceColor = [0.93 0.69 0.13]; b5(3).FaceColor = [0.85 0.33 0.1];
ylabel('Reactions Added'); title('2. Minimal Media Gap (Glucose Stress Test)');
xticklabels(T.Cluster); xtickangle(45); grid on;

nexttile;
bar_data_var = [T.Var_COV, T.Var_FAST, T.Var_iMAT];
b2 = bar(bar_data_var, 'grouped');
b2(1).FaceColor = [0 0.45 0.74]; b2(2).FaceColor = [0.93 0.69 0.13]; b2(3).FaceColor = [0.85 0.33 0.1];
ylabel('Variance (%)'); title('3. Variance Captured by Fixed Functional Models');
xticklabels(T.Cluster); xtickangle(45); grid on; ylim([0 100]);

nexttile;
bar_data_avg = [T.AvgVar_COV, T.AvgVar_FAST, T.AvgVar_iMAT];
b3 = bar(bar_data_avg, 'grouped');
b3(1).FaceColor = [0 0.45 0.74]; b3(2).FaceColor = [0.93 0.69 0.13]; b3(3).FaceColor = [0.85 0.33 0.1];
ylabel('Variance/Rxn'); title('4. Variance Density');
xticklabels(T.Cluster); xtickangle(45); grid on;

nexttile;
bar_data_lost = [T.Lost_COV, T.Lost_FAST, T.Lost_iMAT];
b4 = bar(bar_data_lost, 'grouped');
b4(1).FaceColor = [0 0.45 0.74]; b4(2).FaceColor = [0.93 0.69 0.13]; b4(3).FaceColor = [0.85 0.33 0.1];
ylabel('Lost Reactions'); title('5. Total Lost Reactions per Method');
xticklabels(T.Cluster); xtickangle(45); grid on;

exportgraphics(fig, outputPdfPath);
fprintf('\nReport saved to: %s\n', outputPdfPath);

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function [gap_size, added_indices] = calculate_gap_and_patch(model, active_mask, target_idx)
    [m, n] = size(model.S);
    
    target_flux = 1.0; 
    
    % STRICT MODE: Everything inactive must be paid for.
    penalizable_mask = ~active_mask;
    penalizable_mask(target_idx) = false; 
    
    inactive_idx = find(penalizable_mask);
    k = length(inactive_idx);
    if k == 0, gap_size = 0; added_indices = []; return; end
    
    num_vars = n + k;
    intcon = (n + 1) : num_vars;
    c = [zeros(n, 1); ones(k, 1)]; 
    
    Aeq = [model.S, sparse(m, k)];
    beq = zeros(m, 1);
    
    num_ineq = 1 + 2*k;
    A_ineq = sparse(num_ineq, num_vars);
    b_ineq = zeros(num_ineq, 1);
    
    A_ineq(1, target_idx) = -1;
    b_ineq(1) = -target_flux;
    
    lb_total = [model.lb; zeros(k, 1)];
    ub_total = [model.ub; ones(k, 1)];
    
    for i = 1:k
        rxn_id = inactive_idx(i);
        M_eff = max([1000, abs(model.ub(rxn_id)), abs(model.lb(rxn_id))]); 
        
        A_ineq(1 + i, rxn_id) = 1;
        A_ineq(1 + i, n + i) = -M_eff;
        
        A_ineq(1 + k + i, rxn_id) = -1;
        A_ineq(1 + k + i, n + i) = -M_eff;
    end
    
    opts_milp = optimoptions('intlinprog', 'Display', 'off', ...
        'MaxTime', 120, ... 
        'IntegerTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6);
    
    try
        [x, ~, exitflag] = intlinprog(c, intcon, A_ineq, b_ineq, Aeq, beq, lb_total, ub_total, opts_milp);
        
        if exitflag == 1 || exitflag == 2
            v_vals = x(1:n);
            z_vals = round(x(n+1:end));
            
            % =========================================================
            % THE FIX: Catch the Big-M switch OR any trace leak down to 1e-12!
            % =========================================================
            flux_mask = (z_vals == 1) | (abs(v_vals(inactive_idx)) > 1e-12);
            
            added_indices = inactive_idx(flux_mask);
            gap_size = length(added_indices);
        else
            gap_size = NaN; added_indices = [];
        end
    catch
        gap_size = NaN; added_indices = [];
    end
end
function [kept_indices] = run_imat_milp_standard(model, RH, RL)
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
    if isempty(gpr_indices), RH_indices=[]; RL_indices=[]; return; end
    gpr_expr = rxnExpression(gpr_indices);
    high_thresh = prctile(gpr_expr, 75);
    low_thresh = prctile(gpr_expr, 25);

    fprintf("High: %d, low: %d \n",high_thresh,low_thresh);
    RH_tmp = gpr_indices(gpr_expr >= high_thresh);
    RL_tmp = gpr_indices(gpr_expr <= low_thresh);
    RH_indices = RH_tmp(:); RL_indices = RL_tmp(:);
    fprintf("lenght RH indices: %d \n",numel(RH_indices));
end


function [intact_precursors, dead_precursors] = perform_biomass_autopsy(model, kept_indices, bio_idx)
    % =====================================================================
    % BIOMASS AUTOPSY: Tests which specific precursors can still be 
    % synthesized by the pruned metabolic network.
    % =====================================================================
    
    % 1. Identify Precursors (Negative coefficients in the Biomass reaction)
    bio_col = model.S(:, bio_idx);
    precursor_mask = bio_col < -1e-6; 
    precursor_idx = find(precursor_mask);
    precursor_names = model.mets(precursor_idx);
    num_precursors = length(precursor_idx);
    
    % 2. Apply the Pruning (Shut down the lost reactions)
    [n_mets, n_rxns] = size(model.S);
    pruned_model = model;
    lost_mask = true(n_rxns, 1);
    lost_mask(kept_indices) = false;
    
    pruned_model.lb(lost_mask) = 0;
    pruned_model.ub(lost_mask) = 0;
    
    % 3. Initialize Tracking
    intact_precursors = {};
    dead_precursors = {};
    opts = optimoptions('linprog', 'Display', 'none');
    
    % 4. Loop through every single precursor
    for p = 1:num_precursors
        m_idx = precursor_idx(p);
        m_name = precursor_names{p};
        
        % Create a temporary model with a new Demand reaction for this precursor
        % Equation: 1 * Precursor -> [Out]
        temp_model = pruned_model;
        temp_model.S = [temp_model.S, sparse(n_mets, 1)];
        temp_model.S(m_idx, end) = -1; % Drain the metabolite
        
        % Objective: Maximize this single demand reaction
        temp_model.c = zeros(n_rxns + 1, 1);
        temp_model.c(end) = 1; 
        
        temp_model.lb = [temp_model.lb; 0];
        temp_model.ub = [temp_model.ub; 1000];
        
        % Run FBA
        try
            [~, fval, exitflag] = linprog(-temp_model.c, [], [], temp_model.S, zeros(n_mets,1), temp_model.lb, temp_model.ub, opts);
            
            % If solver succeeds and flux is greater than a tiny threshold, it survives!
            if exitflag == 1 && -fval > 1e-5
                intact_precursors{end+1} = m_name;
            else
                dead_precursors{end+1} = m_name;
            end
        catch
            dead_precursors{end+1} = m_name;
        end
    end
end