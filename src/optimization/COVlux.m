% =========================================================================
% COVLUX: Covariance-based Reconstruction of Elementary Flux Modes
% =========================================================================
% Uses config.json for paths/params and dynamic relative path handling
% =========================================================================

%% 1. SETUP ENVIRONMENT & PATHS
% A. Resolve Paths Relative to This Script
currentScriptPath = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(currentScriptPath)); % src/optimization -> src -> root

% Add Helper Functions
functionsPath = fullfile(currentScriptPath, 'functions');
if ~exist(functionsPath, 'dir')
    error('Functions folder not found at: %s', functionsPath);
end
addpath(functionsPath);

% B. Load Configuration
configFile = fullfile(projectRoot, 'config', 'config.json');
if ~exist(configFile, 'file')
    error('Config file not found at: %s', configFile);
end
config = jsondecode(fileread(configFile));

% --- PARAMETERS FROM CONFIG ---
usebigbasis      = config.params.use_big_basis;
SNR_THRESHOLD    = config.params.snr_threshold;
FLUX_NOISE_FLOOR = config.params.flux_noise_floor;
lambda_l21       = config.params.lambda_l21;
lambda_qr        = config.params.lambda_qr;
max_iters        = config.params.max_iters;
verbose          = config.params.verbose;
div_by_reactions = config.params.div_by_reactions;
mean_influence  = config.params.mean_influence  ;
fprintf("DIV by reactions??: %d \n",div_by_reactions)
second_moment_mat= config.params.second_moment_mat;

% Additional Hardcoded Params
reduce_colinearity    = false;
reduce_EFM_only       = false;
reduce_reactions_only = false;
lambda_lasso          = NaN;
lambda_ridge          = NaN;
lambda_balance        = NaN; 


% --- RESOLVE INPUT PATHS (Relative to Project Root) ---
modelPath   = fullfile(projectRoot, config.paths.models_dir, config.model.model_file);
efmMatPath  = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{1});
efmMatPath2 = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{2});
efmMatPath3 = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{3});
efmMatPath4 = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{4});

% --- DYNAMIC INPUT DIRECTORY (MRAS OUTPUTS) ---
clusteringBase = fullfile(projectRoot, config.paths.results_dir, 'Clustering');
run_name       = config.params.input_clustering_folder;

if isfield(config.params, 'use_conditions') && config.params.use_conditions
    % Pattern: results/Clustering/<Run>/data_files/grouped_by_condition/MRAS_outputs
    covDir = fullfile(clusteringBase, run_name, 'data_files', 'grouped_by_condition', 'MRAS_outputs');
else
    % Pattern: results/Clustering/<Run>/data_files/MRAS_outputs
    covDir = fullfile(clusteringBase, run_name, 'data_files', 'MRAS_outputs');
end

% Fallback: If MRAS_outputs doesn't exist, try the parent folder
if ~exist(covDir, 'dir')
    warning('MRAS_outputs folder not found at: %s\nFalling back one level up.', covDir);
    covDir = fileparts(covDir); 
end
if ~exist(covDir, 'dir')
    error('Input directory not found: %s', covDir);
end
fprintf('COVLUX Input Source: %s\n', covDir);

%% 2. LOAD MODEL & EFM BASIS
fprintf('Loading Metabolic Model...\n');
data = load(modelPath);
if isfield(data, 'pruned_ir')
    model_ir = data.pruned_ir;
else
    vars = fieldnames(data);
    model_ir = data.(vars{1});
end
S_model = model_ir.S;
rxnNames = model_ir.rxns;
[m_met, n_rxn] = size(S_model);

fprintf('Loading EFM Basis Matrices...\n');
S1 = load(efmMatPath, 'EFM_matrix', 'rxnNames');
S2 = load(efmMatPath2, 'EFM_matrix', 'rxnNames');
S4 = load(efmMatPath4, 'EFM_matrix', 'rxnNames');
flux1 = S1.EFM_matrix(2:end, :); res1 = S1.EFM_matrix(1, :);
flux2 = S2.EFM_matrix(2:end, :); res2 = S2.EFM_matrix(1, :);
flux4 = S4.EFM_matrix(2:end, :); res4 = S4.EFM_matrix(1, :);
name1 = string(S1.rxnNames(:));
name2 = string(S2.rxnNames(:));
name4 = string(S4.rxnNames(:));

% Find Intersection of Reactions across bases
common_names = intersect(name1, name2, 'stable');
common_names = intersect(common_names, name4, 'stable');
idx3 = [];
if usebigbasis
    fprintf('Loading Big Basis (S3)...\n');
    S3 = load(efmMatPath3, 'new_EFM_matrix', 'rxnE');
    flux3 = S3.new_EFM_matrix(2:end, :);
    res3 = S3.new_EFM_matrix(1, :);
    name3 = string(S3.rxnE(:));
    
    [rxnE, ~, ~] = intersect(common_names, name3, 'stable');
    [~, idx3] = ismember(rxnE, name3);
else
    rxnE = common_names;
end
[~, idx1] = ismember(rxnE, name1);
[~, idx2] = ismember(rxnE, name2);
[~, idx4] = ismember(rxnE, name4);

% Merge EFM Matrices
fprintf('Merging EFM Bases...\n');
if usebigbasis
    E_combined_raw = [[res1; flux1(idx1,:)], ...
                      [res2; flux2(idx2,:)], ...
                      [res3; flux3(idx3,:)], ...
                      [res4; flux4(idx4,:)]];
else
    E_combined_raw = [[res1; flux1(idx1,:)], ...
                      [res2; flux2(idx2,:)], ...
                      [res4; flux4(idx4,:)]];
end

% Quality Control (Filter Noise)
[E_full_clean, ~] = quality_filter(E_combined_raw, SNR_THRESHOLD, FLUX_NOISE_FLOOR);
E_full = E_full_clean(2:end, :); % Remove resistance row
fprintf('Final Valid EFM Count: %d\n', size(E_full, 2));

rxn_names_full=rxnNames;
% Check Coverage
is_active = any(abs(E_full) > FLUX_NOISE_FLOOR, 2);
active_names = rxnE(is_active);
[~, idx_in_model] = ismember(active_names, rxnNames);
full_coverage = false(n_rxn, 1);
full_coverage(idx_in_model(idx_in_model > 0)) = true;
uncovered_list = rxnNames(~full_coverage);
if ~isempty(uncovered_list)
    fprintf('WARNING: %d Reactions in model are UNCOVERED by EFMs.\n', length(uncovered_list));
end

%% 3. PREPARE OUTPUT DIRECTORIES
if usebigbasis
    baseResultsDir = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_bigbasis');
else
    baseResultsDir = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_smallbasis');
end

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
resultsDir = fullfile(baseResultsDir, ['Run_' timestamp]);
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
plotDir = fullfile(resultsDir, 'plots');
if ~exist(plotDir, 'dir'), mkdir(plotDir); end

copyfile(configFile, fullfile(resultsDir, 'run_settings.json'));

%% 4. PROCESS CLUSTERS (MAIN LOOP)
% Find all COV files recursively or in specific dir
files = dir(fullfile(covDir, '**', 'cluster_*_COV.csv')); 
files = files(~[files.isdir]);
% Dedup logic
[~, uIdx] = unique(lower(fullfile({files.folder}, {files.name})));
files = files(uIdx);

if isempty(files)
    error('No covariance files found in: %s', covDir);
end

allResults = [];
summaryMetrics = [];
for k = 1:numel(files)
    filePath = fullfile(files(k).folder, files(k).name);
    [~, fileName, ~] = fileparts(files(k).name);
    clusterName = regexprep(fileName, '_COV.*$', '');
    
    fprintf('\n=== Processing %d/%d: %s ===\n', k, numel(files), clusterName);
    
    % --- Load Data ---
    try
        T_cov = readtable(filePath, 'TextType', 'string', 'ReadRowNames', true, 'PreserveVariableNames', true);
        % Handle cases where RowNames aren't detected automatically
        if isempty(T_cov.Properties.RowNames)
             T_cov = readtable(filePath, 'TextType', 'string', 'PreserveVariableNames', true);
             rn = string(T_cov{:,1});
             T_cov = T_cov(:, 2:end);
             T_cov.Properties.RowNames = rn;
        end
    catch ME
        warning('Skipping %s due to read error: %s', clusterName, ME.message);
        continue;
    end
    
    outer_path = strrep(filePath, '_COV.csv', '_MEAN_OUTER.csv');
    if exist(outer_path, 'file')
        T_outer = readtable(outer_path, 'TextType', 'string', 'ReadRowNames', true, 'PreserveVariableNames', true);
    else
        T_outer = [];
    end
    

    second_moment_path = strrep(filePath, '_COV.csv', '_SecondMoment.csv');


    if second_moment_mat
        if exist(second_moment_path, 'file')
            % 1. Read the Table
            T_temp = readtable(second_moment_path, 'TextType', 'string', 'ReadRowNames', true, 'PreserveVariableNames', true);
            
            % 2. Convert to Numeric Matrix IMMEDIATELY
            X_full = table2array(T_temp);
            
            % 3. (Optional) Extract Reaction Names from the table rows if needed
            rxnX = string(T_temp.Properties.RowNames);
        else
            [X_full, rxnX, ~, lambda_balance] = align_and_combine_C_and_M(T_cov, T_outer, plotDir, clusterName,verbose);
        end
    else 
        [X_full, rxnX, ~, lambda_balance] = align_and_combine_C_and_M(T_cov, T_outer, plotDir, clusterName,verbose);
    end


    % Intersect with EFM Matrix
    [rxnC, iE, iX] = intersect(rxnE, rxnX, 'stable');
    
    if isempty(rxnC)
        warning('No overlap between EFMs and Data for %s', clusterName);
        continue;
    end
    E_original_full=E_full;
    E_use = E_full(iE, :);
    X_use = X_full(iX, iX);
    rxns_current = rxnE(iE);
    
    % --- Sanitize Covariance ---
    % Enforce PSD and clean noise
    [X_use, stats] = sanitize_covariance_matrix(X_use, 1e-5, verbose, plotDir, clusterName);
    %[X_use, scale_factor] = prepare_covariance_matrix(X_use,1000,verbose);
    
    % --- Load Mean Activity (mu) ---
    % Attempt to find corresponding MEAN file
    mean_filename = [clusterName '_MEAN.csv'];
    % Look in same folder or parent
    mean_path = fullfile(files(k).folder, mean_filename);
    if ~exist(mean_path, 'file')
        % Try removing '_MRAS' or other suffixes if naming varies
        mean_path = fullfile(files(k).folder, strrep(mean_filename, '_MRAS', ''));
    end
    
    mu_use = zeros(length(rxns_current), 1);
    if exist(mean_path, 'file')
        T_mean = readtable(mean_path, 'ReadRowNames', true, 'VariableNamingRule', 'preserve');
        % Assume values are in first column
        mean_vals = table2array(T_mean(:,1));
        mean_map = containers.Map(string(T_mean.Properties.RowNames), mean_vals);
        
        for r = 1:length(rxns_current)
            if isKey(mean_map, rxns_current(r))
                mu_use(r) = mean_map(rxns_current(r));
            end
        end
    else
        % Fallback: approximation from diagonal
        mu_use = sqrt(diag(X_use)); 
    end
    
    % --- Optimization Prep ---
    E_final = E_use;
    X_final = X_use;
    mu_final = mu_use;
    
    P_full = eye(size(E_use, 1));
    P_efm = eye(size(E_use, 2));
    
    % Calculate 'Reduced Mean' (Projection of mean onto EFM space)
    mu_final_reduced = abs(E_final' * mu_final);
    
    % --- RUN OPTIMIZATION (Local Function) ---
    protected_idx = [];
    %[A_opt_QR, E_red, L] = covlux_symmetric_pipeline_mean(E_final, X_final, mu_final, lambda_qr, lambda_l21, max_iters, mean_influence, protected_idx, plotDir, clusterName,verbose,div_by_reactions);
    %[A_opt_QR, E_red, L] = covlux_symmetric_pipeline_mean_lasso(E_final, X_final, mu_final, lambda_qr, lambda_l21, max_iters, mean_influence, protected_idx, plotDir, clusterName,verbose,div_by_reactions);
    %[A_opt_QR, E_red, L, metrics]= solve_weighted_lasso_covariance(E_final, X_final, verbose, plotDir, clusterName);
    [A_opt_QR, E_red, L, metrics]= covariance_selection(E_final, X_final, verbose, plotDir, clusterName,mean_influence,lambda_l21);
    tolerance = 1e-9; 
    
    efm_norms_local = sqrt(sum(E_final.^2, 1));
    efm_norms_local(efm_norms_local < 1e-12) = 1; % Avoid div-by-zero
    E_norm_check = E_final ./ efm_norms_local;

    
    X_recon_full = E_norm_check * A_opt_QR * E_norm_check';

    
    recon_error = norm(X_final - X_recon_full, 'fro') / norm(X_final, 'fro');
    
    fprintf('  Reconstruction Error: %.6f (%.2f%%)\n', recon_error, recon_error * 100);
    
    fprintf("Smallest value in A: %f",min(A_opt_QR(:)));
    A_full = A_opt_QR;
    
    variances = diag(A_full);
    selected_efm_idx = find(variances > 1e-9);
   
   
    X_recon_reduced = X_recon_full;  
    
    cluster_metrics = compute_and_save_cluster_metrics(...
        resultsDir, clusterName, ...
        E_use, X_use, E_final, X_final, A_opt_QR, P_full, P_efm, ...
        reduce_colinearity, reduce_EFM_only, reduce_reactions_only, ...
        lambda_lasso, lambda_ridge, X_recon_full, L, E_red, selected_efm_idx, rxns_current, ...
        E_original_full, rxn_names_full, ...
        lambda_l21, max_iters, mean_influence, mu_final,lambda_balance);

    
    % Add to aggregate results
    if isempty(allResults)
        allResults = cluster_metrics;
    else
        allResults(end+1) = cluster_metrics;
    end
    
    if isempty(summaryMetrics)
        summaryMetrics = cluster_metrics;
    else
        summaryMetrics(end+1) = cluster_metrics;
    end
    
end

if exist('allResults', 'var') && ~isempty(allResults)
    results_to_save = allResults;
elseif exist('summaryMetrics', 'var') && ~isempty(summaryMetrics)
    results_to_save = summaryMetrics;
else
    fprintf('No results to save.\n');
    return;
end

% Convert to table
final_summary_table = struct2table(results_to_save);

% Save as single CSV
final_summary_file = fullfile(resultsDir, 'all_clusters_summary.csv');
writetable(final_summary_table, final_summary_file);

% PRINT THE ENTIRE TABLE
fprintf('\n=== COMPREHENSIVE RESULTS SUMMARY ===\n');
disp(final_summary_table);

fprintf('\n✓ Saved comprehensive summary: %s\n', final_summary_file);
fprintf('  - %d clusters\n', height(final_summary_table));
fprintf('  - %d metrics per cluster\n', width(final_summary_table));



%% ========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [A_final, E_reduced, L_final] = covlux_symmetric_pipeline_mean(E, X, mu_efm, lambda_qr, lambda_l21, max_iters, mean_influence, protected_idx, saveDir, clusterName,verbose,div_by_reactions)
    % COVLUX with Rent-Based Selection & Full Visualization Suite
    
    
    [m, s] = size(E);
    if verbose, fprintf('  [Normalization] Scaling all EFMs to Unit Norm to prevent magnitude bias...\n'); end
    
    E_orig=E;
    
    efm_norms = sqrt(sum(E.^2, 1));
    
    
    efm_norms(efm_norms < 1e-12) = 1; 
    
    
    E = E ./ efm_norms;
    % ============================================================
    % PLOTTING BLOCK 1: Scale Analysis
    % ============================================================
    fprintf('Saving optimisation figures to: %s\n', saveDir);

    fig1 = figure('Name', 'Scale Mismatch Analysis', 'Color', 'w', 'Visible', 'off');
    
    subplot(2, 1, 1);
    histogram(E(:), 100, 'FaceColor', [0 0.4470 0.7410]);
    title('Distribution of E (Normalized Input)');
    xlabel('Value'); ylabel('Count (Log Scale)');
    set(gca, 'YScale', 'log'); grid on;
    
    subplot(2, 1, 2);
    histogram(X(:), 100, 'FaceColor', [0.8500 0.3250 0.0980]);
    title('Distribution of X_{orig} (Target Data)');
    xlabel('Value'); ylabel('Count (Log Scale)');
    set(gca, 'YScale', 'log'); grid on;
    
    if nargin >= 9 && ~isempty(saveDir)
        filename = fullfile(saveDir, sprintf('%s_scale_mismatch.png', clusterName));
        saveas(fig1, filename);
    end
    close(fig1);
    
    % ============================================================
    % INITIAL SOLVE
    % ============================================================
   
    A_opt = solve_QR_regularized_new(E, X, lambda_qr);
    A = A_opt;
    
    
    E_reaction_norms = sum(E.^2, 1)'; 
    
    % ============================================================
    % PLOTTING BLOCK 2: Initial Distribution
    % ============================================================
    row_strengths = sqrt(sum(A.^2, 2)); 
    
    h_diag = figure('Name', 'Initial Dense A Distribution', 'Color', 'w', 'Visible', 'off');
    
    subplot(2, 1, 1);
    histogram(row_strengths, 200, 'FaceColor', [0.8500 0.3250 0.0980]);
    set(gca, 'YScale', 'log'); 
    title(sprintf('Distribution of EFM Strengths (Max=%.2e)', max(row_strengths)));
    xlabel('Row Strength (L2 Norm)'); ylabel('Count (Log Scale)'); grid on;
    
    sorted_str = sort(row_strengths, 'descend');
    subplot(2, 1, 2);
    plot(sorted_str, 'LineWidth', 3, 'Color', 'k');
    title('Scree Plot: EFM Strength Ranked High to Low');
    xlabel('Rank'); ylabel('Strength');
    xlim([1, min(length(sorted_str), 2000)]); grid on;
    
    if nargin >= 9 && ~isempty(saveDir)
        filename = fullfile(saveDir, sprintf('%s_initial_distribution.png', clusterName));
        saveas(h_diag, filename);
    end
    close(h_diag);
    
    % ============================================================
    % OPTIMIZATION: ADAPTIVE BACKWARD SELECTION
    % ============================================================
    fprintf('\n--- Starting Adaptive Backward Selection ---\n');
    
    % Metric Storage
    qr_error_history = nan(max_iters, 1); 
    sparsity_history = zeros(max_iters, 1);
    iteration_numbers = zeros(max_iters, 1);
    
    history_k = [];
    history_r2 = [];
    history_score = [];
    history_A = {}; 
    
    active_mask = true(s, 1);
    current_k = s;
    P_candidates = size(E, 2); 

    

    for iter = 1:max_iters
        %current_variances = diag(A) .* E_reaction_norms;

        current_variances = sum(A,2) .* E_reaction_norms;
        if div_by_reactions
            
            reactions_per_efm = sum(abs(E) > 1e-6, 1)';
            ranking_scores = current_variances ./ (reactions_per_efm + 1e-6);
        else 
            ranking_scores = current_variances;
        end
        ranking_scores(~active_mask) = inf; 
        if ~isempty(protected_idx)
            ranking_scores(protected_idx) = inf;
        end
        
        [sorted_scores, sort_idx] = sort(ranking_scores, 'ascend');
        candidates = sort_idx(isfinite(sorted_scores));
        
        if isempty(candidates), break; end
        
        % Adaptive Batch Size
        max_batch_size = ceil(current_k * 0.05); 
        max_batch_size = max(1, max_batch_size); 
        n_safe = min(length(candidates), max_batch_size);
        potential_victims = candidates(1:n_safe);
        victim_scores = sorted_scores(1:n_safe);
        
        % Cliff Detection
        n_to_kill = 1; 
        if length(potential_victims) > 1
            jump_threshold = 5.0; 
            for k = 2:length(potential_victims)
                curr_score = victim_scores(k);
                prev_score = victim_scores(k-1);
                if prev_score < 1e-15
                    if curr_score > 1e-15, break; else, n_to_kill = k; end
                else
                    if (curr_score / prev_score) > jump_threshold, break; else, n_to_kill = k; end
                end
            end
        end
        
        kill_idx = candidates(1:n_to_kill);
        active_mask(kill_idx) = false;
        current_k = sum(active_mask);
        
        % Refit
        if current_k > 0
            E_reduced_iter = E(:, active_mask);
            A_small = solve_QR_regularized_new(E_reduced_iter, X, lambda_qr);
            A_next = zeros(s, s);
            idx_active = find(active_mask);
            A_next(idx_active, idx_active) = A_small;
            A = A_next;
        else
            A = zeros(s, s);
        end
        
        % Metrics
        X_rec = E * A * E';
        qr_fit_error = norm(X_rec - X, 'fro') / norm(X, 'fro');
        current_r2 = max(0, 1 - qr_fit_error^2);
        current_sparsity = 1 - (current_k / s);
        
        
        if current_k > 0
            active_efms = find(active_mask);
        
           
            reactions_used = any(abs(E(:, active_efms)) > 0, 2);
            
            n_active_reactions = sum(reactions_used);
        else
            n_active_reactions = 0;
        end
        
        
        m_total = size(X,1); 
        
        
        k_params = current_k; 
        n_reactions_lost=m_total - n_active_reactions;
        
        reaction_sparsity_score = n_reactions_lost / m_total;

        %efficiency_score = (1 - lambda_l21) * current_r2 + lambda_l21 * reaction_sparsity_score;

        % Or multiplication version (more aggressive trade-off)
        efficiency_score = current_r2 * (reaction_sparsity_score ^ lambda_l21);

        % Store History
        qr_error_history(iter) = qr_fit_error;
        sparsity_history(iter) = current_sparsity;
        iteration_numbers(iter) = iter;
        
        history_k = [history_k; current_k]; 
        history_r2 = [history_r2; current_r2];
        history_score = [history_score; efficiency_score];
        history_A{end+1} = A; 
        
        if verbose
            fprintf('  Iter %d: Killed %d | R^2=%f | k=%d | active Rxns=%d | Score=%.6e\n', ...
                iter, n_to_kill, current_r2 ,current_k, n_active_reactions, efficiency_score);
        else
            if mod(iter, 10) == 0
                fprintf('  Iter %d: Killed %d | R^2=%f | k=%d | active Rxns=%d | Score=%.6e\n', ...
                    iter, n_to_kill, current_r2 ,current_k, n_active_reactions, efficiency_score);
            end 
        end 
        if current_r2 < 0.01 || current_k == 0
            break;
        end
    end
    
    % ============================================================
    % SELECTION: THE "QUALITY FLOOR" FIX
    % ============================================================
    
    
    
    min_r2_threshold = 0.1; 
    valid_mask = history_r2 > min_r2_threshold;
    
    if ~any(valid_mask)
        fprintf('  WARNING: No models met the R2 > %.2f threshold. Picking best available.\n', min_r2_threshold);
        [~, best_idx] = max(history_score);
    else
        
        valid_indices = find(valid_mask);
        valid_scores = history_score(valid_indices);
        [~, local_max_idx] = max(valid_scores);
        
        best_idx = valid_indices(local_max_idx);
    end
    
    best_iter = best_idx;
    best_score = history_score(best_idx);
    best_A = history_A{best_idx};
    
    fprintf('  Optimal Solution Selected: Iter %d (k=%d, R2=%.4f, Score=%.6e)\n', ...
        best_iter, history_k(best_iter), history_r2(best_iter), best_score);
    % ============================================================
    % PLOTTING BLOCK 3: Convergence History (Restored)
    % ============================================================
    fig2 = figure('Name', 'Convergence History', 'Color', 'w', 'Visible', 'off');
    
    subplot(2, 1, 1);
    valid_idx = iteration_numbers > 0;
    x_qr = iteration_numbers(valid_idx);
    y_qr = qr_error_history(valid_idx);
    
    plot(x_qr, y_qr, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'g');
    xlabel('Iteration'); ylabel('QR Refit Error');
    title('Convergence: Error vs Iteration');
    grid on; ylim([0, 1]); 
    
    subplot(2, 1, 2);
    plot(iteration_numbers(valid_idx), sparsity_history(valid_idx), 'r-', 'LineWidth', 2);
    xlabel('Iteration'); ylabel('Sparsity');
    title('Sparsity Evolution');
    grid on; ylim([0, 1]);
    
    if nargin >= 9 && ~isempty(saveDir)
        filename = fullfile(saveDir, sprintf('%s_convergence_history.png', clusterName));
        saveas(fig2, filename);
        
    end
    close(fig2);

    % ============================================================
    % PLOTTING BLOCK 4: The Scree / Elbow Plot (New)
    % ============================================================
    fig_elbow = figure('Name', 'Optimization Scree Plot', 'Color', 'w', 'Visible', 'off');
    
    % Top Panel: R2 vs k (The Elbow)
    subplot(2, 1, 1);
    plot(history_k, history_r2, '-bo', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    hold on;
    opt_k = history_k(best_iter);
    opt_r2 = history_r2(best_iter);
    plot(opt_k, opt_r2, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); 
    text(opt_k, opt_r2 - 0.05, sprintf(' Optimal Rent (k=%d)', opt_k), 'Color', 'r', 'FontWeight', 'bold');
    hold off;
    
    set(gca, 'XDir', 'reverse'); 
    xlabel('Number of Active EFMs (k)'); 
    ylabel('Explained Variance (R^2)');
    title('Scree Plot: R^2 vs Model Complexity');
    grid on;
    
    % Bottom Panel: Efficiency Score
    subplot(2, 1, 2);
    plot(history_k, history_score, '-k', 'LineWidth', 2);
    hold on;
    xline(opt_k, '--r', 'Max Efficiency');
    hold off;
    set(gca, 'XDir', 'reverse');
    xlabel('Number of Active EFMs (k)');
    ylabel('Efficiency Score (R^2 / \surd{k})');
    title('Optimization Function (Rent Sweet Spot)');
    grid on;
    
    if nargin >= 9 && ~isempty(saveDir)
        filename = fullfile(saveDir, sprintf('%s_scree_optimization.png', clusterName));
        saveas(fig_elbow, filename);
        
    end
    close(fig_elbow);
    
    % ============================================================
    % FINAL SELECTION & RECONSTRUCTION
    % ============================================================
    
    % Use the 'best_A' determined by the Rent Score
    final_row_norms = sqrt(sum(best_A.^2, 2));
    selected_idx = find(final_row_norms > 1e-9);
    
    if ~isempty(protected_idx)
        selected_idx = unique([selected_idx; protected_idx]);
        selected_idx = sort(selected_idx); 
    end
    
    fprintf('  Final Selection: %d of %d EFMs (%.1f%%)\n', length(selected_idx), s, 100*length(selected_idx)/s);
    
    E_reduced = E_orig(:, selected_idx);
    
    % Plot length distribution of survivors
    plot_efm_length_distribution(E_reduced, saveDir, clusterName);
    
    % Final Polish
    [A_reduced, L_final] = solve_psd_final_qr_psd(E_reduced, X,verbose,saveDir,clusterName);
    A_final = zeros(s, s);
    A_final(selected_idx, selected_idx) = A_reduced;
    
    % ============================================================
    % PLOTTING BLOCK 5: Final Solution Analysis
    % ============================================================
    fig4 = figure('Name', 'Final Solution Analysis', 'Color', 'w', 'Visible', 'off');
    
    subplot(2, 2, 1);
    A_nonzero = A_final(A_final ~= 0);
    if ~isempty(A_nonzero)
        histogram(log10(abs(A_nonzero)), 50, 'FaceColor', [0.2 0.6 0.2]);
        xlabel('log_{10}(|A_{ij}|)'); ylabel('Count'); title('Non-zero entries of A'); grid on;
    end
    
    subplot(2, 2, 2);
    final_row_norms_A = sqrt(sum(A_final.^2, 2));
    histogram(log10(final_row_norms_A(final_row_norms_A > 0)), 30, 'FaceColor', [0.8 0.2 0.2]);
    xlabel('log_{10}(Row Norm)'); ylabel('Count'); title('EFM Row Norms'); grid on;
    
    subplot(2, 2, 3);
    bar(sort(final_row_norms_A, 'descend'));
    xlabel('EFM Index (sorted)'); ylabel('Row Norm'); title('EFM Importance');
    set(gca, 'YScale', 'log'); grid on;
    
    subplot(2, 2, 4);
    selected_counts = sum(abs(E(:, selected_idx)) > 0, 2);
    bar(sort(selected_counts, 'descend'));
    xlabel('Reaction Index (sorted)'); ylabel('Participating EFMs'); title('Reaction Participation'); grid on;
    
    if nargin >= 9 && ~isempty(saveDir)
        filename = fullfile(saveDir, sprintf('%s_final_solution_analysis.png', clusterName));
        saveas(fig4, filename);
        
    end
    close(fig4);
    
    final_error = norm(E * A_final * E' - X, 'fro') / norm(X, 'fro');
    fprintf("final error: %.f")

end



function [A_opt_QR, E_red, L] = covlux_symmetric_pipeline_mean_lasso(E, X, mu_efm, lambda_qr, lambda_l21, max_iters, mean_influence, protected_idx, plotDir, clusterName, verbose, div_by_reactions)


    [m, s] = size(E);
    
    
    if verbose, fprintf('  [Normalization] Scaling all EFMs to Unit Rays...\n'); end
    
    E_orig = E; 
    
    efm_norms = sqrt(sum(E.^2, 1));
   
    efm_norms(efm_norms < 1e-12) = 1; 
    
    
    E_norm = E ./ efm_norms; 
    
    [A_final, E_reduced_raw, L_final_raw] = solve_weighted_lasso_covariance(E_orig, X, verbose, plotDir, clusterName);
    
    A_opt_QR = A_final;
   
    selected_idx = find(diag(A_opt_QR) > 1e-12);
    
    E_red = E_norm(:, selected_idx);
 
    active_activities = sqrt(diag(A_opt_QR(selected_idx, selected_idx)));
    L = active_activities; 
    
    X_recon = E_norm * A_opt_QR * E_norm';
   
    recon_error = norm(X - X_recon, 'fro') / norm(X, 'fro');
    
    fprintf('  Final Reconstruction Error: %.4f (%.2f%%)\n', recon_error, recon_error*100);
    
    
    mag_ratio = sum(X_recon(:)) / sum(X(:));
    fprintf('  Magnitude Check (Recon/Data): %.4f (Should be close to 1.0)\n', mag_ratio);
end