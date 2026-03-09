% =========================================================================
% BIOLOGICAL STATE ANALYSIS (Batch Mode for Latest Run)
% =========================================================================
% 1. Dynamically finds the project root.
% 2. Auto-detects the LATEST Run directory.
% 3. Loads config from the run's archived run_settings.json.
% 4. Dynamically merges and filters raw EFM bases exactly like COVLUX.
% 5. Calculates Phenotypes, Objectives, Bottlenecks, and Redundancy.
% 6. Visualizes Pruning Capacity Trajectories (Subsystems, Energy, AA, MTA)
% =========================================================================
TARGET_K = []; % Leave empty [] to use the dynamic Optimal K per cluster
%% 1. ROBUST PROJECT ROOT FINDER & AUTO-DETECT LATEST RUN
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
    error('Could not locate "config/config.json". Please run from project root.');
end
% Add required functions to path (for quality_filter)
functionsPath = fullfile(projectRoot, 'src', 'optimization', 'functions');
if exist(functionsPath, 'dir')
    addpath(functionsPath);
else
    warning('Functions path not found at: %s. quality_filter may fail.', functionsPath);
end
% We use the master config strictly to locate the base results folder
masterConfigFile = fullfile(projectRoot, 'config', 'config.json');
masterConfig = jsondecode(fileread(masterConfigFile));
if masterConfig.params.use_big_basis
    covluxBase = fullfile(projectRoot, masterConfig.paths.results_dir, 'COVlux_cov_bigbasis');
else
    covluxBase = fullfile(projectRoot, masterConfig.paths.results_dir, 'COVlux_cov_smallbasis');
end
d = dir(fullfile(covluxBase, 'Run_*'));
d = d([d.isdir]);
if isempty(d)
    error('No COVLUX run folders found in: %s', covluxBase);
end
[~, idx] = sort([d.datenum], 'descend');
latestRun = d(idx(1)).name;
runDir = fullfile(covluxBase, latestRun);
fprintf('\nTargeting LATEST COVLUX Results: %s\n', runDir);
logFileName = fullfile(runDir, sprintf('Analysis_Log_%s.txt', datestr(now, 'yyyy-mm-dd_HH-MM')));
diary(logFileName); 
fprintf('\n=======================================================\n');
fprintf('LOGGING STARTED: %s\n', datestr(now));
fprintf('LOG FILE: %s\n', logFileName);
fprintf('Targeting LATEST COVLUX Results: %s\n', runDir);
fprintf('=======================================================\n');
%% 2. LOAD STRICTLY FROM ARCHIVED RUN SETTINGS
settingsFile = fullfile(runDir, 'run_settings.json');
if ~exist(settingsFile, 'file')
    error('run_settings.json not found in the run directory: %s', runDir);
end
config = jsondecode(fileread(settingsFile)); 
run_name         = config.params.input_clustering_folder;
lambda_l21       = config.params.lambda_l21;
mean_influence   = config.params.mean_influence;
usebigbasis      = config.params.use_big_basis;
SNR_THRESHOLD    = config.params.snr_threshold;
FLUX_NOISE_FLOOR = config.params.flux_noise_floor;
fprintf('\n=======================================================\n');
fprintf('                MODEL RUN PARAMETERS                   \n');
fprintf('=======================================================\n');
fprintf('  Project Name (Cells): %s\n', run_name);
fprintf('  Lambda L2,1:          %.4f\n', lambda_l21);
fprintf('  Mean Influence:       %.4f\n', mean_influence);
if isempty(TARGET_K)
    fprintf('  Evaluation K:         [Dynamic] Optimal K per cluster\n');
else
    fprintf('  Evaluation K:         %d EFMs (Hardcoded)\n', TARGET_K);
end
fprintf('=======================================================\n\n');
%% 3. LOAD MODEL & DYNAMICALLY BUILD EFM BASIS
% --- RESOLVE INPUT PATHS (Relative to Project Root) ---
modelPath   = fullfile(projectRoot, config.paths.models_dir, config.model.model_file);
efmMatPath  = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{1});
efmMatPath2 = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{2});
efmMatPath3 = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{3});
efmMatPath4 = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{4});
fprintf('Loading Metabolic Model: %s\n', modelPath);
data = load(modelPath);
if isfield(data, 'pruned_ir')
    model_ir = data.pruned_ir;
else
    vars = fieldnames(data);
    model_ir = data.(vars{1});
end
S_model = model_ir.S;
rxnNames = model_ir.rxns;
model = model_ir; % Map for downstream subsystem lookup
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
n_total_efms = size(E_full, 2);
fprintf('Final Valid EFM Count: %d\n', n_total_efms);
rxn_names_full = rxnNames;
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
%% 4. FIND ALL CLUSTERS AND LOOP
log_dir = fullfile(runDir, 'log_dir');
if ~exist(log_dir, 'dir'), mkdir(log_dir); end
logFiles = dir(fullfile(log_dir, '*_removal_log.csv'));
if isempty(logFiles)
    error('No removal logs (*_removal_log.csv) found in %s', log_dir);
end
fprintf('\nFound %d cluster(s) to analyze.\n', length(logFiles));
for f = 1:length(logFiles)
    clusterName = strrep(logFiles(f).name, '_logCPM_removal_log.csv', '');
    clusterName = strrep(clusterName, '_logCPM', '');
    logPath = fullfile(log_dir, logFiles(f).name);
    
    fprintf('\n\n#######################################################\n');
    fprintf('              ANALYZING: %s\n', upper(clusterName));
    fprintf('#######################################################\n');
    
    % Read the specific log
    removal_log = readtable(logPath);
    
    % --- DETERMINE EVALUATION K ---
    if isempty(TARGET_K)
        current_k = removal_log.Optimal_K(1);
        fprintf('Targeting Optimal Solver Cutoff: K = %d\n', current_k);
    else
        current_k = TARGET_K;
        fprintf('Targeting Hardcoded Cutoff: K = %d\n', current_k);
    end
    
    % --- RECONSTRUCT SURVIVING EFMs ---
    cluster_initial_ids = unique(removal_log.Removed_Global_ID);
    removed_before_target = removal_log.Removed_Global_ID(removal_log.Remaining_K >= current_k);
    surviving_ids = setdiff(cluster_initial_ids, removed_before_target);
    E_survivors = E_full(:, surviving_ids);
    
    fprintf('=== RECONSTRUCTED STATE AT K = %d ===\n', length(surviving_ids));
    
    % --- 1. I/O Phenotype (Exchange Reactions) ---
    fprintf('\n--- 1. I/O Phenotype (Exchange Reactions) ---\n');
    ex_idx = find(startsWith(string(rxnE), 'EX_', 'IgnoreCase', true));
    
    if ~isempty(ex_idx)
        ex_usage = sum(abs(E_survivors(ex_idx, :)) > 1e-6, 2) / length(surviving_ids);
        [sorted_usage, sort_order] = sort(ex_usage, 'descend');
        top_ex_idx = ex_idx(sort_order);
        
        fprintf('Top active exchanges (Present in %% of active EFMs):\n');
        for i = 1:min(10, length(top_ex_idx))
            if sorted_usage(i) > 0
                fprintf('  %s: %.1f%%\n', rxnE(top_ex_idx(i)), sorted_usage(i)*100);
            end
        end
    end
    
    % --- 2. Functional Objectives ---
    fprintf('\n--- 2. Functional Objectives ---\n');
    obj_keywords = {'biomass', 'growth', 'ATPM', 'maintenance'};
    for k = 1:length(obj_keywords)
        idx = find(contains(string(rxnE), obj_keywords{k}, 'IgnoreCase', true));
        for i = 1:length(idx)
            usage = sum(abs(E_survivors(idx(i), :)) > 1e-6) / length(surviving_ids);
            fprintf('  %s: %.1f%% of active EFMs\n', rxnE(idx(i)), usage*100);
        end
    end
    
    % --- 3. Bottleneck Analysis ---
    fprintf('\n--- 3. Core Bottlenecks & Essentiality ---\n');
    rxn_usage = sum(abs(E_survivors) > 1e-6, 2) / length(surviving_ids);
    internal_idx = setdiff(1:length(rxnE), ex_idx);
    internal_usage = rxn_usage(internal_idx);
    
    [sorted_int_usage, int_sort_order] = sort(internal_usage, 'descend');
    top_int_idx = internal_idx(int_sort_order);
    
    fprintf('Top Internal Bottlenecks (Strictly required by this state):\n');
    count = 0;
    for i = 1:length(top_int_idx)
        if sorted_int_usage(i) >= 0.90
            count = count + 1;
            rxn_name = rxnE(top_int_idx(i));
            
            % Look up the reaction in the original model to get its Subsystem mapping
            mod_idx = find(strcmp(model.rxns, rxn_name));
            pathway_str = '';
            if ~isempty(mod_idx) && isfield(model, 'subSystems') && ~isempty(model.subSystems{mod_idx})
                pathway = model.subSystems{mod_idx};
                if iscell(pathway), pathway = pathway{1}; end
                pathway_str = sprintf(' [%s]', char(pathway));
            end
            
            fprintf('  %s%s: %.1f%%\n', rxn_name, pathway_str, sorted_int_usage(i)*100);
        end
        if count >= 15, break; end 
    end
    
    % --- 4. Final State Overlap & Redundancy ---
    fprintf('\n--- 4. Final State Overlap & Redundancy ---\n');
    binary_E = abs(E_survivors) > 1e-6;
    if size(binary_E, 2) > 1
        jaccard_dist = pdist(binary_E', 'jaccard');
        avg_jaccard = mean(1 - jaccard_dist);
        
        cosine_dist = pdist(E_survivors', 'cosine');
        avg_cosine = mean(1 - cosine_dist);
        
        fprintf('  Average Pairwise Jaccard Similarity: %.2f%%\n', avg_jaccard * 100);
        fprintf('  Average Pairwise Cosine Similarity:  %.2f%%\n', avg_cosine * 100);
        
        if avg_jaccard > 0.85
            fprintf('  -> INSIGHT: High redundancy. Surviving EFMs are highly similar variations of the same core pathway.\n');
        elseif avg_jaccard < 0.40
            fprintf('  -> INSIGHT: High diversity. The cell relies on structurally distinct, parallel metabolic strategies.\n');
        else
            fprintf('  -> INSIGHT: Balanced state with a mix of shared core reactions and unique peripheral routes.\n');
        end
    else
        fprintf('  Only 1 EFM left. No pairwise comparison possible.\n');
    end
    % =====================================================================
    % 5. TRAJECTORY VISUALIZATION (PRUNING CAPACITY)
    % =====================================================================
    fprintf('\n--- Generating Capacity Trajectory Visualizations ---\n');
    
    % Define the timeline (X-axis: K values from max down to target)
    max_k = max(removal_log.Remaining_K);
    k_steps = unique(round(linspace(max_k, current_k, 50)), 'stable');
    k_steps = sort(k_steps, 'descend'); % Time flows from left (max_k) to right (current_k)
    
    % --- Mapping Subsystems ---
    subsys_map = containers.Map();
    for i = 1:length(rxnE)
        mod_idx = find(strcmp(model.rxns, rxnE(i)));
        if ~isempty(mod_idx) && isfield(model, 'subSystems') && ~isempty(model.subSystems{mod_idx})
            pw = model.subSystems{mod_idx};
            if iscell(pw), pw = pw{1}; end
            if ~isempty(pw)
                if ~isKey(subsys_map, char(pw))
                    subsys_map(char(pw)) = [];
                end
                subsys_map(char(pw)) = [subsys_map(char(pw)), i];
            end
        end
    end
    subsys_names = keys(subsys_map);
    
    % --- Category Definitions for Heatmaps ---
    % Define keyword sets for each category
    energy_keywords = {'glc','o2','atp','adp','nad','nadh','fad','fdh','atpm'};
    aa_keywords     = {'trp','met','cys','glu','gln','ala','arg','asn','asp','his','ile','leu','lys','phe','pro','ser','thr','tyr','val','nh4'};
    mta_keywords    = {'mta','mtan','sah','sam','methan','methyl'};
    % I/O: any reaction starting with 'EX_' (exchange)
    
    % Initialize category index lists
    energy_idx = [];
    aa_idx     = [];
    mta_idx    = [];
    io_idx     = [];
    
    % Scan all reactions and assign to categories
    for i = 1:length(rxnE)
        rxn_str = lower(string(rxnE(i)));
        % I/O (exchange) reactions
        if startsWith(rxn_str, 'ex_')
            io_idx(end+1) = i;
        else
            % Check other categories
            if any(contains(rxn_str, energy_keywords))
                energy_idx(end+1) = i;
            end
            if any(contains(rxn_str, aa_keywords))
                aa_idx(end+1) = i;
            end
            if any(contains(rxn_str, mta_keywords))
                mta_idx(end+1) = i;
            end
        end
    end
    
    
    % Create a cell array of category data for easy looping
    categories = {
        'Energy', energy_idx;
        'Amino_Acids', aa_idx;
        'MTA_Cofactors', mta_idx;
        'IO_Exchanges', io_idx
    };
    
    % --- Calculate Capacities over Time ---
    E_init_bin = abs(E_full(:, cluster_initial_ids)) > 1e-6;
    
    % Baseline capacity (Denominator)
    subsys_base = zeros(length(subsys_names), 1);
    for s = 1:length(subsys_names)
        subsys_base(s) = sum(E_init_bin(subsys_map(subsys_names{s}), :), 'all');
    end
    
    % Track capacities: subsystems and categories
    subsys_cap_traj = zeros(length(subsys_names), length(k_steps));
    % For categories, we'll store a structure with matrix of size (n_reactions_in_category x time)
    cat_traj = cell(size(categories,1),1);
    cat_rxns = cell(size(categories,1),1);
    
    for t = 1:length(k_steps)
        step_k = k_steps(t);
        step_removed = removal_log.Removed_Global_ID(removal_log.Remaining_K > step_k);
        step_survivors = setdiff(cluster_initial_ids, step_removed);
        
        if isempty(step_survivors)
            continue; 
        end
        
        E_step_bin = abs(E_full(:, step_survivors)) > 1e-6;
        
        % Subsystems
        for s = 1:length(subsys_names)
            curr_val = sum(E_step_bin(subsys_map(subsys_names{s}), :), 'all');
            if subsys_base(s) > 0
                subsys_cap_traj(s, t) = (curr_val / subsys_base(s)) * 100;
            end
        end
        
        % Categories
        for c = 1:size(categories,1)
            idx_list = categories{c,2};
            if isempty(idx_list), continue; end
            % For each reaction in this category, compute remaining capacity
            % We'll store a matrix of size (length(idx_list) x length(k_steps))
            if t == 1
                % Pre-allocate
                cat_traj{c} = zeros(length(idx_list), length(k_steps));
                cat_rxns{c} = rxnE(idx_list);
                % Baseline for normalization
                cat_base = sum(E_init_bin(idx_list, :), 2); % per reaction total counts
                cat_base(cat_base == 0) = 1; % avoid division by zero
                cat_traj{c}(:,t) = (sum(E_step_bin(idx_list, :), 2) ./ cat_base) * 100;
            else
                cat_base = sum(E_init_bin(idx_list, :), 2);
                cat_base(cat_base == 0) = 1;
                cat_traj{c}(:,t) = (sum(E_step_bin(idx_list, :), 2) ./ cat_base) * 100;
            end
        end
    end
    
    
    % --- PLOT 1: Subsystem Shrinkage (Top 10 pruned) with Tiled Layout ---
    % Find subsystems that lost the most capacity
    cap_loss = 100 - subsys_cap_traj(:, end);
    [~, sort_loss] = sort(cap_loss, 'descend');
    top_plot_idx = sort_loss(1:min(10, length(sort_loss)));
    
    f1 = figure('Name', 'Subsystem Capacity Pruning', 'Position', [100, 100, 1200, 500], 'Visible', 'off');
    tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Left tile: Line chart
    nexttile;
    hold on; grid on;
    colors = lines(length(top_plot_idx));
    for i = 1:length(top_plot_idx)
        plot(k_steps, subsys_cap_traj(top_plot_idx(i), :), 'LineWidth', 2.5, 'Color', colors(i,:), ...
            'DisplayName', subsys_names{top_plot_idx(i)});
    end
    set(gca, 'XDir', 'reverse');
    xlabel('Number of Surviving EFMs (K)');
    ylabel('Remaining Capacity (%)');
    title(sprintf('Pathway Pruning (%s)', upper(clusterName)));
    legend('Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 8);
    
    % Right tile: Pie chart of remaining capacity at final K
    nexttile;
    final_cap = subsys_cap_traj(:, end);
    % Exclude subsystems with zero final capacity or zero base? but final_cap could be zero.
    % Take top 5 subsystems by remaining capacity, group rest.
    [sorted_cap, sort_idx] = sort(final_cap, 'descend');
    top_n = 30;
    if length(sorted_cap) > top_n
        top_vals = sorted_cap(1:top_n);
        top_names = subsys_names(sort_idx(1:top_n));
        other_val = sum(sorted_cap(top_n+1:end));
        if other_val > 0
            top_vals(end+1) = other_val;
            top_names{end+1} = 'Others';
        end
    else
        top_vals = sorted_cap;
        top_names = subsys_names(sort_idx);
    end
    % Only plot if there are positive values
    if any(top_vals > 0)
        pie(top_vals, top_names);
        title(sprintf('Final Capacity Distribution (K=%d)', current_k));
    else
        text(0.5,0.5,'No remaining capacity','HorizontalAlignment','center');
        title('Final Capacity');
    end
    plotDir = fullfile(runDir, 'plots'); 
    if ~exist(plotDir, 'dir')
        mkdir(plotDir);
    end
    saveas(f1, fullfile(plotDir, sprintf('%s_Subsystem_Trajectory.png', clusterName)));
    close(f1);

    % --- PLOT 1B: DECOMPRESSED PIE CHART (INDEPENDENT LINES) ---
    top_n_lines = 30; % Keeping exactly the 30 categories you requested!
    
    top_idx_pie = sort_idx(1:min(top_n_lines, length(sort_idx)));
    other_idx_pie = sort_idx(top_n_lines+1:end);
    
    pie_traj_data = zeros(length(top_idx_pie) + 1, length(k_steps));
    pie_traj_names = cell(length(top_idx_pie) + 1, 1);
    
    % Populate data exactly from the capacities used in the Pie Chart
    for i = 1:length(top_idx_pie)
        pie_traj_data(i, :) = subsys_cap_traj(top_idx_pie(i), :);
        pie_traj_names{i} = strrep(subsys_names{top_idx_pie(i)}, '_', ' ');
    end
    
    % Sum the remaining minor subsystems into a single "Others" block
    if ~isempty(other_idx_pie)
        pie_traj_data(end, :) = sum(subsys_cap_traj(other_idx_pie, :), 1);
    else
        pie_traj_data(end, :) = 0;
    end
    pie_traj_names{end} = 'Others';
    
    % Normalize to 100% at every K step (Relative Share of the Pie)
    pie_traj_perc = (pie_traj_data ./ sum(pie_traj_data, 1)) * 100;
    pie_traj_perc(isnan(pie_traj_perc)) = 0; % Prevent NaN errors
    
    % Create a wider figure to accommodate a 30-item legend
    f1b = figure('Name', 'Subsystem Share Over Time', 'Position', [150, 150, 1600, 700], 'Visible', 'off');
    hold on; grid on;
    
    % Use a colormap designed to generate distinct colors for large arrays
    line_colors = turbo(length(top_idx_pie)); 
    line_colors = [line_colors; 0.7 0.7 0.7]; % Add gray for the "Others" line
    
    % Plot INDEPENDENT LINES instead of stacked areas so they don't push each other up!
    for i = 1:length(pie_traj_names)
        plot(k_steps, pie_traj_perc(i, :), 'LineWidth', 2.5, 'Color', line_colors(i,:), 'DisplayName', pie_traj_names{i});
    end
    
    set(gca, 'XDir', 'reverse'); 
    
    xlabel('Number of Surviving EFMs (K)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Relative Share of the Network (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Independent Subsystem Shares Over Time (Decompressed) (%s)', upper(clusterName)), 'FontSize', 14);
    
    % Set legend to 2 columns so all 30 items fit nicely on the screen
    legend('Location', 'eastoutside', 'FontSize', 9, 'Interpreter', 'none', 'NumColumns', 2);
    
    saveas(f1b, fullfile(plotDir, sprintf('%s_1B_Subsystem_Distribution_Lines.png', clusterName)));
    close(f1b);
    % --- SURVIVORS VS. LOSERS ---
    
    
    % 1. Filter out subsystems that had zero baseline capacity
    valid_sys_mask = subsys_base > 0;
    valid_indices = find(valid_sys_mask);
    valid_final_cap = subsys_cap_traj(valid_sys_mask, end);
    
    % 2. Sort by how much capacity they retained at the end
    [sorted_cap, sort_idx] = sort(valid_final_cap, 'descend');
    
    % 3. Define how many to compare (e.g., Top 10 of each)
    num_compare = min(10, floor(length(sort_idx)/2)); 
    
    % Grab the extreme ends of the sorted list
    survivor_idx_local = sort_idx(1:num_compare);        % Highest final capacity
    loser_idx_local    = sort_idx(end-num_compare+1:end); % Lowest final capacity
    
    survivor_idx_global = valid_indices(survivor_idx_local);
    loser_idx_global    = valid_indices(loser_idx_local);
    
    % 4. Create the Figure
    f_compare = figure('Name', 'Survivors vs Losers', 'Position', [200, 200, 1100, 600], 'Visible', 'off');
    hold on; grid on;
    
    % 5. Plot the SURVIVORS (Solid, thick, cool colors)
    survivor_colors = winter(num_compare); % Blues/Greens
    for i = 1:num_compare
        idx = survivor_idx_global(i);
        plot(k_steps, subsys_cap_traj(idx, :), 'LineWidth', 3, 'Color', survivor_colors(i,:), ...
            'DisplayName', sprintf('[KEPT] %s (%.1f%%)', subsys_names{idx}, subsys_cap_traj(idx, end)));
    end
    
    % 6. Plot the LOSERS (Dashed, thinner, warm colors)
    loser_colors = autumn(num_compare); % Reds/Oranges/Yellows
    for i = 1:num_compare
        % We iterate backwards so the absolute worst loser gets plotted first in deep red
        idx = loser_idx_global(num_compare - i + 1); 
        plot(k_steps, subsys_cap_traj(idx, :), 'LineWidth', 2, 'LineStyle', '--', 'Color', loser_colors(i,:), ...
            'DisplayName', sprintf('[DROPPED] %s (%.1f%%)', subsys_names{idx}, subsys_cap_traj(idx, end)));
    end
    
    % 7. Formatting
    set(gca, 'XDir', 'reverse'); % K goes from max down to target
    xlabel('Number of Surviving EFMs (K)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Remaining Capacity (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('The Core Shift: Highly Preserved vs. Abandoned Pathways (%s)', upper(clusterName)), 'FontSize', 14);
    
    % Add a nice legend outside the plot area
    legend('Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 10);
    
    % Save it
    saveas(f_compare, fullfile(plotDir, sprintf('%s_Survivors_vs_Losers.png', clusterName)));
    close(f_compare);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % =====================================================================
    % --- NEW ADDITION: RELATIVE ENRICHMENT & DUAL ANALYSIS ---
    % =====================================================================
    fprintf('  Calculating Relative Enrichment metrics...\n');
    
    subsys_cap_traj_rel = zeros(length(subsys_names), length(k_steps)); 
    
    for t = 1:length(k_steps)
        step_k = k_steps(t);
        step_removed = removal_log.Removed_Global_ID(removal_log.Remaining_K > step_k);
        step_survivors = setdiff(cluster_initial_ids, step_removed);
        
        if isempty(step_survivors), continue; end
        
        E_step_bin = abs(E_full(:, step_survivors)) > 1e-6;
        
        for s = 1:length(subsys_names)
            rxn_idx = subsys_map(subsys_names{s});
            % RELATIVE ENRICHMENT (Funnel)
            efms_using_sys = any(E_step_bin(rxn_idx, :), 1);
            if step_k > 0
                subsys_cap_traj_rel(s, t) = (sum(efms_using_sys) / step_k) * 100;
            end
        end
    end
    
    % --- DUAL PLOT: SURVIVORS VS. LOSERS (Waterfall vs Funnel) ---
    fprintf('  Generating Dual Comparison Plot (Pruning vs Enrichment)...\n');
    
    valid_sys_mask = subsys_base > 0;
    valid_indices = find(valid_sys_mask);
    valid_final_rel = subsys_cap_traj_rel(valid_sys_mask, end);
    
    [~, sort_idx_rel] = sort(valid_final_rel, 'descend');
    num_compare_dual = min(7, floor(length(sort_idx_rel)/2)); 
    
    surv_idx_dual  = valid_indices(sort_idx_rel(1:num_compare_dual));
    loser_idx_dual = valid_indices(sort_idx_rel(end-num_compare_dual+1:end));
    
    f_dual = figure('Name', 'Dual Trajectory Analysis', 'Position', [100, 100, 1600, 700], 'Visible', 'off');
    t_layout = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    surv_colors_dual = winter(num_compare_dual); 
    loser_colors_dual = autumn(num_compare_dual);
    
    % LEFT TILE: Absolute Pruning (The Waterfall)
    ax1 = nexttile; hold on; grid on;
    for i = 1:num_compare_dual
        idx = surv_idx_dual(i);
        plot(k_steps, subsys_cap_traj(idx, :), 'LineWidth', 3, 'Color', surv_colors_dual(i,:));
    end
    for i = 1:num_compare_dual
        idx = loser_idx_dual(num_compare_dual - i + 1); 
        plot(k_steps, subsys_cap_traj(idx, :), 'LineWidth', 2, 'LineStyle', '--', 'Color', loser_colors_dual(i,:));
    end
    set(ax1, 'XDir', 'reverse');
    xlabel(ax1, 'Number of Surviving EFMs (K)', 'FontWeight', 'bold');
    ylabel(ax1, 'Absolute Remaining Capacity (%)', 'FontWeight', 'bold');
    title(ax1, 'Physical Pruning', 'FontSize', 14);
    
    % RIGHT TILE: Relative Enrichment (The Funnel)
    ax2 = nexttile; hold on; grid on;
    for i = 1:num_compare_dual
        idx = surv_idx_dual(i);
        plot(k_steps, subsys_cap_traj_rel(idx, :), 'LineWidth', 3, 'Color', surv_colors_dual(i,:), ...
            'DisplayName', sprintf('[KEPT] %s (%.1f%%)', subsys_names{idx}, subsys_cap_traj_rel(idx, end)));
    end
    for i = 1:num_compare_dual
        idx = loser_idx_dual(num_compare_dual - i + 1); 
        plot(k_steps, subsys_cap_traj_rel(idx, :), 'LineWidth', 2, 'LineStyle', '--', 'Color', loser_colors_dual(i,:), ...
            'DisplayName', sprintf('[DROPPED] %s (%.1f%%)', subsys_names{idx}, subsys_cap_traj_rel(idx, end)));
    end
    set(ax2, 'XDir', 'reverse');
    xlabel(ax2, 'Number of Surviving EFMs (K)', 'FontWeight', 'bold');
    ylabel(ax2, 'Prevalence in Surviving Pool (%)', 'FontWeight', 'bold');
    title(ax2, 'Biological Enrichment ', 'FontSize', 14);
    
    legend(ax2, 'Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 9);
    title(t_layout, sprintf('The Core Shift: Highly Preserved vs. Abandoned Pathways (%s)', upper(clusterName)), 'FontSize', 16, 'FontWeight', 'bold');
    
    saveas(f_dual, fullfile(plotDir, sprintf('%s_Dual_Trajectory_Analysis.png', clusterName)));
    close(f_dual);
    % --- NEW PLOT: THE BIGGEST MOVERS (Largest Change in Enrichment) ---
    fprintf('  Generating "Biggest Movers" Enrichment Plot...\n');
    
    % Calculate the absolute change in relative enrichment (Final - Initial)
    enrichment_change = subsys_cap_traj_rel(:, end) - subsys_cap_traj_rel(:, 1);
    
    % Find the biggest growers (positive change) and biggest crashers (negative change)
    [~, sort_change] = sort(enrichment_change, 'descend');
    num_movers = 5; % Top 5 up, Top 5 down
    
    biggest_gainers = sort_change(1:num_movers);
    biggest_losers  = sort_change(end-num_movers+1:end);
    
    f_movers = figure('Name', 'Biggest Movers', 'Position', [150, 150, 1000, 600], 'Visible', 'off');
    hold on; grid on;
    
    % Plot Gainers (Thick Green/Blue)
    gainer_colors = winter(num_movers);
    for i = 1:num_movers
        idx = biggest_gainers(i);
        plot(k_steps, subsys_cap_traj_rel(idx, :), 'LineWidth', 3, 'Color', gainer_colors(i,:), ...
            'DisplayName', sprintf('[SURGED +%.1f%%] %s', enrichment_change(idx), subsys_names{idx}));
    end
    
    % Plot Losers (Thick Red/Orange)
    loser_colors = autumn(num_movers);
    for i = 1:num_movers
        idx = biggest_losers(num_movers - i + 1); % Reverse to get the worst at the bottom
        plot(k_steps, subsys_cap_traj_rel(idx, :), 'LineWidth', 3, 'LineStyle', '-.', 'Color', loser_colors(i,:), ...
            'DisplayName', sprintf('[CRASHED %.1f%%] %s', enrichment_change(idx), subsys_names{idx}));
    end
    
    set(gca, 'XDir', 'reverse');
    xlabel('Number of Surviving EFMs (K)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Relative Enrichment (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('The Biggest Movers: Highest Change in Relative Enrichment (%s)', upper(clusterName)), 'FontSize', 14);
    legend('Location', 'eastoutside', 'Interpreter', 'none', 'FontSize', 10);
    
    saveas(f_movers, fullfile(plotDir, sprintf('%s_Biggest_Movers_Enrichment.png', clusterName)));
    close(f_movers);
    % =====================================================================
    
    % --- PLOT 2-5: Category Heatmaps ---
    for c = 1:size(categories,1)
        cat_name = categories{c,1};
        traj_data = cat_traj{c};
        rxn_list = cat_rxns{c};
        
        if isempty(traj_data) || size(traj_data,1) == 0
            fprintf('  Skipping %s heatmap (no reactions in this category).\n', cat_name);
            continue;
        end
        
        f2 = figure('Name', sprintf('%s Heatmap', cat_name), 'Position', [150, 150, 1000, max(400, 20*size(traj_data,1))], 'Visible', 'off');
        imagesc(traj_data);
        colormap(parula);
        c = colorbar;
        c.Label.String = 'Capacity Remaining (%)';
        
        % Formatting
        set(gca, 'YTick', 1:length(rxn_list), 'YTickLabel', rxn_list, 'TickLabelInterpreter', 'none', 'FontSize', 8);
        
        % Simplify X-axis ticks (Show K values)
        num_ticks = min(10, length(k_steps));
        tick_idx = round(linspace(1, length(k_steps), num_ticks));
        set(gca, 'XTick', tick_idx, 'XTickLabel', k_steps(tick_idx));
        
        xlabel('Number of Surviving EFMs (K)');
        title(sprintf('%s Pruning (%s)', strrep(cat_name,'_',' '), upper(clusterName)));
        
        saveas(f2, fullfile(plotDir, sprintf('%s_%s_Heatmap.png', clusterName, cat_name)));
        close(f2);
    end
    
    fprintf('  -> Saved capacity trajectory plots to: %s\n', plotDir);
    % =====================================================================
    % 6. FBA MAX PRODUCTION SIMULATION (TRUE BIOLOGICAL MINIMAL CAPACITY)
    % =====================================================================
    fprintf('\n--- Simulating FBA Max Production on Reduced Networks ---\n');
    
    max_k_val = max(removal_log.Remaining_K);
    k_timeline = unique(round(linspace(max_k_val, current_k, 50)), 'stable');
    k_timeline = sort(k_timeline, 'descend'); 
    % 1. Setup Base Bounds
    if isfield(model, 'lb'), lb_base = model.lb; else, lb_base = zeros(n_rxn, 1); end
    if isfield(model, 'ub'), ub_base = model.ub; else, ub_base = 1000 * ones(n_rxn, 1); end
    if isfield(model, 'rev')
        for idx_r=1:n_rxn, if model.rev(idx_r)==0 && lb_base(idx_r)<0, lb_base(idx_r)=0; end, end
    end
    % 2. ENFORCE TRUE BIOLOGICAL MINIMAL MEDIUM
    % Close ALL uptake reactions first (starve the cell completely)
    uptake_idx = find(startsWith(lower(rxnNames), 'ex_') & endsWith(lower(rxnNames), '_b'));
    lb_base(uptake_idx) = 0; 
    ub_base(uptake_idx) = 0; 
    
    % Open ONLY the scientifically required minimal substrates for E. coli
    minimal_substrates = {'ex_glc__d_e_b', 'ex_glc_e_b', 'ex_o2_e_b', 'ex_h2o_e_b', 'ex_h_e_b', ...
                          'ex_nh4_e_b', 'ex_pi_e_b', 'ex_so4_e_b', 'ex_k_e_b', 'ex_na1_e_b', ...
                          'ex_mg2_e_b', 'ex_ca2_e_b', 'ex_cl_e_b', 'ex_fe2_e_b', 'ex_fe3_e_b', ...
                          'ex_zn2_e_b', 'ex_mn2_e_b', 'ex_cu2_e_b', 'ex_cobalt2_e_b', ...
                          'ex_mobd_e_b', 'ex_ni2_e_b', 'ex_sel_e_b', 'ex_tungs_e_b', ...
                          'ex_thm_e_b', 'ex_cbl1_e_b', 'ex_niacin_e_b'};
                          
    for i = 1:length(minimal_substrates)
        idx = find(strcmpi(rxnNames, minimal_substrates{i}));
        if ~isempty(idx)
            ub_base(idx) = 1000; % Open the door for essential trace elements
        end
    end
    
    % Put realistic physiological limits on Carbon and Oxygen to prevent math exploits
    glc_idx = find(strcmpi(rxnNames, 'ex_glc__d_e_b') | strcmpi(rxnNames, 'ex_glc_e_b'));
    if ~isempty(glc_idx), ub_base(glc_idx) = 10; end 
    o2_idx = find(strcmpi(rxnNames, 'ex_o2_e_b'));
    if ~isempty(o2_idx), ub_base(o2_idx) = 20; end
    % 3. Define Targets: Demand (DM), Export (EX_..._f), Biomass, ATPM
    target_rxns_mod = {};
    for i = 1:length(rxnNames)
        r_name = lower(rxnNames{i});
        is_output = startsWith(r_name, 'ex_') || startsWith(r_name, 'dm_') || contains(r_name, 'biomass') || contains(r_name, 'atpm');
        has_kw = any(contains(r_name, {'trp','met','cys','glu','gln','ala','arg','asn','asp','his','ile','leu','lys','phe','pro','ser','thr','tyr','val','mta', 'atp', 'nadh', 'fad'}));
        if is_output && has_kw && ~endsWith(r_name, '_b')
            target_rxns_mod{end+1} = rxnNames{i};
        end
    end
    target_rxns_mod = unique(target_rxns_mod)';
    
    if ~isempty(target_rxns_mod)
        prod_traj = zeros(length(target_rxns_mod), length(k_timeline));
        options = optimoptions('linprog', 'Display', 'off'); 
        Aeq = model.S; beq = zeros(m_met, 1);
        
        for t = 1:length(k_timeline)
            step_ids = setdiff(cluster_initial_ids, removal_log.Removed_Global_ID(removal_log.Remaining_K > k_timeline(t)));
            if isempty(step_ids), continue; end
            
            % Generate topological mask based strictly on EFM active state
            active_in_E = any(abs(E_full(:, step_ids)) > 1e-6, 2);
            [~, active_mod_idx] = ismember(rxnE(active_in_E), rxnNames);
            active_mod_idx = active_mod_idx(active_mod_idx > 0);
            
            inactive_mask = true(n_rxn, 1); inactive_mask(active_mod_idx) = false;
            
            for p = 1:length(target_rxns_mod)
                t_idx = find(strcmp(rxnNames, target_rxns_mod{p}), 1);
                
                lb_step = lb_base; ub_step = ub_base;
                
                % Apply EFM topological knockout
                lb_step(inactive_mask) = 0; ub_step(inactive_mask) = 0;
                
                % Protect the specific exhaust pipe so we can measure output
                lb_step(t_idx) = 0; ub_step(t_idx) = 1000; 
                
                f = zeros(n_rxn, 1); f(t_idx) = -1; % Maximize output
                [~, fval, exitflag] = linprog(f, [], [], Aeq, beq, lb_step, ub_step, options);
                
                if exitflag == 1
                    prod_traj(p, t) = -fval; 
                else
                    prod_traj(p, t) = 0; 
                end
            end
        end
        
        % 4. Normalize to Initial 100% Capacity
        prod_traj_perc = zeros(size(prod_traj));
        for p = 1:length(target_rxns_mod)
            base_prod = prod_traj(p, 1); 
            if base_prod > 1e-6
                prod_traj_perc(p, :) = (prod_traj(p, :) / base_prod) * 100; 
            elseif max(prod_traj(p, :)) > 1e-6
                prod_traj_perc(p, :) = (prod_traj(p, :) / max(prod_traj(p, :))) * 100;
            end
        end
        
        % Only plot targets that had at least *some* internal capacity > 0
        active_prod_mask = max(prod_traj_perc, [], 2) > 0.1;
        if any(active_prod_mask)
            f_prod = figure('Name', 'Physical Production Capacity', 'Position', [150, 150, 1000, max(600, 20*sum(active_prod_mask))], 'Visible', 'off');
            imagesc(prod_traj_perc(active_prod_mask, :)); colormap(parula); cbar = colorbar; cbar.Label.String = 'Relative Internal Max Production (%)';
            set(gca, 'YTick', 1:sum(active_prod_mask), 'YTickLabel', target_rxns_mod(active_prod_mask), 'TickLabelInterpreter', 'none', 'FontSize', 8);
            num_ticks = min(10, length(k_timeline)); tick_idx = round(linspace(1, length(k_timeline), num_ticks));
            set(gca, 'XTick', tick_idx, 'XTickLabel', k_timeline(tick_idx)); xlabel('Surviving EFMs (K)');
            title(sprintf('Minimal Medium Capacity (%s)', upper(clusterName)));
            saveas(f_prod, fullfile(plotDir, sprintf('%s_FBA_Minimal_Capacity.png', clusterName))); close(f_prod);
        else
            fprintf('  -> No internal production capacity found on minimal media. Cell requires rich media supplements.\n');
        end
    end
    % =====================================================================
    % 5. TRAJECTORY VISUALIZATION (MACRO-STATES, MICRO-MOVERS, & DIET)
    % =====================================================================
    fprintf('\n--- Generating Capacity Trajectory Visualizations ---\n');
    
    % Define the timeline (X-axis: K values from max down to target)
    max_k = max(removal_log.Remaining_K);
    k_steps = unique(round(linspace(max_k, current_k, 50)), 'stable');
    k_steps = sort(k_steps, 'descend'); % Time flows from left to right
    
    plotDir = fullfile(runDir, 'plots'); 
    if ~exist(plotDir, 'dir'), mkdir(plotDir); end
    
    % *** MASSIVE SPEEDUP: Precompute the global binary matrix ONCE ***
    E_full_bin = abs(E_full) > 1e-6; 
    
    % --- 5A. GRANULAR MAPPING (For Micro-Movers) ---
    [~, mod_idx_all] = ismember(rxnE, model.rxns);
    subsys_map = containers.Map();
    for i = 1:length(rxnE)
        mod_idx = mod_idx_all(i);
        if mod_idx > 0 && isfield(model, 'subSystems') && ~isempty(model.subSystems{mod_idx})
            pw = model.subSystems{mod_idx};
            if iscell(pw), pw = pw{1}; end
            if ~isempty(pw)
                if ~isKey(subsys_map, char(pw)), subsys_map(char(pw)) = []; end
                subsys_map(char(pw)) = [subsys_map(char(pw)), i];
            end
        end
    end
    subsys_names = keys(subsys_map);

    % --- 5B. EXCLUSIVE MACRO-CATEGORIES (NO DUPLICATES) ---
    % --- 5B. EXCLUSIVE MESO-CATEGORIES (FINER DISTINCTIONS) ---
    fprintf('  Assigning EFMs to exclusive Meso-Categories (Finer Distinctions)...\n');
    super_cat_names = {'Central Carbon', 'Amino Acids', 'Nucleotides', 'Lipids & Envelope', 'Cofactors & Vitamins', 'Transport & Exchange', 'Other'};
    
    kw_central = {'glycolysis', 'gluconeogenesis', 'tca', 'citric', 'pentose phosphate', 'pyruvate', 'carbon', 'glyoxylate'};
    kw_aa      = {'amino acid', 'valine', 'leucine', 'isoleucine', 'alanine', 'aspartate', 'glutamate', 'histidine', 'phenylalanine', 'tyrosine', 'tryptophan', 'arginine', 'proline', 'cysteine', 'methionine', 'threonine', 'serine', 'glycine'};
    kw_nuc     = {'nucleotide', 'purine', 'pyrimidine', 'salvage'};
    kw_lipid   = {'lipid', 'fatty acid', 'envelope', 'murein', 'lps', 'peptidoglycan', 'membrane'};
    kw_cofac   = {'cofactor', 'vitamin', 'folate', 'porphyrin', 'mta', 'methyl', 'nad', 'flavin', 'coenzyme'};
    kw_trans   = {'transport', 'exchange', 'porin', 'efflux', 'extracellular'};

    % Map every reaction to exactly ONE meso-category
    rxn_to_cat = 7 * ones(length(rxnE), 1); % Default is 7 ('Other')
    for i = 1:length(rxnE)
        mod_idx = mod_idx_all(i);
        if mod_idx > 0 && isfield(model, 'subSystems') && ~isempty(model.subSystems{mod_idx})
            pw = lower(char(model.subSystems{mod_idx}));
            if any(contains(pw, kw_central)), rxn_to_cat(i) = 1;
            elseif any(contains(pw, kw_aa)), rxn_to_cat(i) = 2;
            elseif any(contains(pw, kw_nuc)), rxn_to_cat(i) = 3;
            elseif any(contains(pw, kw_lipid)), rxn_to_cat(i) = 4;
            elseif any(contains(pw, kw_cofac)), rxn_to_cat(i) = 5;
            elseif any(contains(pw, kw_trans)), rxn_to_cat(i) = 6;
            end
        end
    end
    
    % DYNAMIC UBIQUITY FILTER
    valid_cats = []; total_efms = size(E_full_bin, 2);
    fprintf('  [Ubiquity Filter Results]:\n');
    for c = 1:7
        in_cat_mask = (rxn_to_cat == c);
        if sum(in_cat_mask) > 0
            prevalence = sum(any(E_full_bin(in_cat_mask, :), 1)) / total_efms;
            if prevalence > 0.98
                fprintf('    -> IGNORING: "%s" (Present in %.1f%% of all EFMs)\n', super_cat_names{c}, prevalence*100);
            else
                fprintf('    -> KEEPING:  "%s" (Present in %.1f%% of all EFMs)\n', super_cat_names{c}, prevalence*100);
                valid_cats(end+1) = c;
            end
        end
    end
    
    % WINNER-TAKES-ALL ASSIGNMENT
    C_mat = zeros(length(valid_cats), length(rxnE));
    for i = 1:length(valid_cats), C_mat(i, rxn_to_cat == valid_cats(i)) = 1; end
    EFM_scores = C_mat * E_full_bin; 
    [max_scores, best_idx] = max(EFM_scores, [], 1);
    
    global_efm_assignments = zeros(1, total_efms);
    for j = 1:total_efms
        if max_scores(j) > 0, global_efm_assignments(j) = valid_cats(best_idx(j)); end
    end
    
    plot_cats = [valid_cats, 0]; 
    plot_cat_names = super_cat_names(valid_cats); 
    plot_cat_names{end+1} = 'Core/Ubiquitous Only';

    % --- 5C. METABOLIC DIET PREP ---
    rxnE_lower = lower(string(rxnE));
    carbon_sources = {'Glucose', 'EX_glc__D_e_b'; 'Acetate', 'EX_ac_e_b'; 'Glycerol', 'EX_glyc_e_b'; 'Succinate', 'EX_succ_e_b'; 'Lactate', 'EX_lac__D_e_b'; 'Fructose', 'EX_fru_e_b'};
    o2_idx = find(rxnE_lower == 'ex_o2_e_b' | rxnE_lower == 'ex_o2_e', 1);
    src_idx = zeros(size(carbon_sources, 1), 1);
    for i = 1:size(carbon_sources, 1)
        base_rxn = lower(carbon_sources{i, 2});
        idx = find(rxnE_lower == base_rxn | rxnE_lower == strrep(base_rxn, '_b', ''), 1);
        if ~isempty(idx), src_idx(i) = idx; end
    end
    active_src_idx = src_idx(src_idx > 0); active_sources = carbon_sources(src_idx > 0, :);

    % --- Pre-allocate Data Arrays ---
    super_counts_abs = zeros(length(plot_cats), length(k_steps)); 
    super_counts_rel = zeros(length(plot_cats), length(k_steps)); 
    subsys_cap_traj_rel = zeros(length(subsys_names), length(k_steps)); 
    src_traj_rel = zeros(length(active_src_idx), length(k_steps));
    strat_aerobic = zeros(1, length(k_steps)); strat_ferment = zeros(1, length(k_steps));

    initial_assignments = global_efm_assignments(cluster_initial_ids);
    initial_totals = zeros(length(plot_cats),1);
    for i = 1:length(plot_cats), initial_totals(i) = sum(initial_assignments == plot_cats(i)); end

    % --- 5D. THE UNIFIED FAST LOOP ---
    fprintf('  Calculating Trajectories...\n');
    for t = 1:length(k_steps)
        step_k = k_steps(t);
        step_removed = removal_log.Removed_Global_ID(removal_log.Remaining_K > step_k);
        step_survivors = setdiff(cluster_initial_ids, step_removed);
        if isempty(step_survivors), continue; end
        
        E_step_bin = E_full_bin(:, step_survivors);
        num_surv = length(step_survivors);
        
        % 1. Macro-Categories (Exclusive)
        surv_assignments = global_efm_assignments(step_survivors);
        for i = 1:length(plot_cats)
            c = plot_cats(i);
            current_count = sum(surv_assignments == c);
            if initial_totals(i) > 0, super_counts_abs(i, t) = (current_count / initial_totals(i)) * 100; end
            if num_surv > 0, super_counts_rel(i, t) = (current_count / num_surv) * 100; end
        end

        % 2. Micro-Movers (Granular Subsystems)
        for s = 1:length(subsys_names)
            efms_using_sys = any(E_step_bin(subsys_map(subsys_names{s}), :), 1);
            if num_surv > 0, subsys_cap_traj_rel(s, t) = (sum(efms_using_sys) / num_surv) * 100; end
        end

        % 3. Metabolic Diet & Strategy
        if num_surv > 0 && (~isempty(active_src_idx) || ~isempty(o2_idx))
            for s = 1:length(active_src_idx)
                src_traj_rel(s, t) = (sum(E_step_bin(active_src_idx(s), :)) / num_surv) * 100;
            end
            
            if ~isempty(active_src_idx), is_carbon_active = any(E_step_bin(active_src_idx, :), 1); else, is_carbon_active = true(1, num_surv); end
            if ~isempty(o2_idx), is_aerobic_efm = E_step_bin(o2_idx, :); else, is_aerobic_efm = false(1, num_surv); end
            
            strat_aerobic(t) = (sum(is_carbon_active & is_aerobic_efm) / num_surv) * 100;
            strat_ferment(t) = (sum(is_carbon_active & ~is_aerobic_efm) / num_surv) * 100;
        end
    end

    % =====================================================================
    % --- PLOT 1: MACRO-BIOLOGY (100% STACKED AREA) ---
    f_macro = figure('Name', 'Macro-State Trajectory', 'Position', [100, 100, 1500, 600], 'Visible', 'off');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Expanded Color Palette for the 7 Categories
    master_colors = containers.Map({1, 2, 3, 4, 5, 6, 7, 0}, ...
        {[0.85, 0.32, 0.09], ... % Red (Central Carbon)
         [0.00, 0.44, 0.74], ... % Blue (Amino Acids)
         [0.30, 0.75, 0.93], ... % Cyan (Nucleotides)
         [0.49, 0.18, 0.56], ... % Purple (Lipids)
         [0.93, 0.69, 0.13], ... % Gold (Cofactors)
         [0.46, 0.67, 0.18], ... % Green (Transport)
         [0.50, 0.50, 0.50], ... % Gray (Other)
         [0.20, 0.20, 0.20]});   % Dark Gray (Core/Ubiquitous)
         
    dyn_colors = zeros(length(plot_cats), 3); 
    for i=1:length(plot_cats), dyn_colors(i,:) = master_colors(plot_cats(i)); end

    % LEFT: Absolute Pruning
    ax1 = nexttile; hold on; grid on;
    for i = 1:length(plot_cats), plot(k_steps, super_counts_abs(i, :), 'LineWidth', 3, 'Color', dyn_colors(i,:)); end
    set(ax1, 'XDir', 'reverse'); xlabel(ax1, 'Surviving EFMs (K)', 'FontWeight', 'bold'); ylabel(ax1, 'Absolute Remaining Count (%)', 'FontWeight', 'bold'); title(ax1, 'Physical EFM Pruning', 'FontSize', 14);

    % RIGHT: Relative Enrichment (100% Stacked Area)
    ax2 = nexttile; hold on; grid on;
    h_area = area(k_steps, super_counts_rel', 'LineStyle', 'none');
    for i = 1:length(plot_cats)
        h_area(i).FaceColor = dyn_colors(i,:); h_area(i).FaceAlpha = 0.8;
        h_area(i).DisplayName = sprintf('%s (%.1f%%)', plot_cat_names{i}, super_counts_rel(i, end));
    end
    set(ax2, 'XDir', 'reverse'); ylim(ax2, [0, 100]);
    xlabel(ax2, 'Surviving EFMs (K)', 'FontWeight', 'bold'); ylabel(ax2, 'Composition of Surviving Pool (%)', 'FontWeight', 'bold'); 
    title(ax2, 'Biological Focus Shift (Exclusive Assignment)', 'FontSize', 14); legend(ax2, 'Location', 'eastoutside', 'FontSize', 10);
    saveas(f_macro, fullfile(plotDir, sprintf('%s_Macro_Trajectory_Analysis.png', clusterName))); close(f_macro);

    % =====================================================================
    % --- PLOT 2: BIGGEST MICRO-MOVERS ---
    enrichment_change = subsys_cap_traj_rel(:, end) - subsys_cap_traj_rel(:, 1);
    [~, sort_change] = sort(enrichment_change, 'descend'); num_movers = 5; 
    
    f_movers = figure('Name', 'Biggest Movers', 'Position', [150, 150, 1000, 600], 'Visible', 'off'); hold on; grid on;
    gainer_colors = winter(num_movers); loser_colors = autumn(num_movers);
    for i = 1:num_movers
        idx = sort_change(i);
        plot(k_steps, subsys_cap_traj_rel(idx, :), 'LineWidth', 3, 'Color', gainer_colors(i,:), 'DisplayName', sprintf('[SURGED +%.1f%%] %s', enrichment_change(idx), subsys_names{idx}));
    end
    for i = 1:num_movers
        idx = sort_change(end-num_movers+i); 
        plot(k_steps, subsys_cap_traj_rel(idx, :), 'LineWidth', 3, 'LineStyle', '-.', 'Color', loser_colors(i,:), 'DisplayName', sprintf('[CRASHED %.1f%%] %s', enrichment_change(idx), subsys_names{idx}));
    end
    set(gca, 'XDir', 'reverse'); xlabel('Surviving EFMs (K)', 'FontSize', 12, 'FontWeight', 'bold'); ylabel('Relative Enrichment (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('The Biggest Micro-Movers (%s)', upper(clusterName)), 'FontSize', 14); legend('Location', 'eastoutside', 'Interpreter', 'none');
    saveas(f_movers, fullfile(plotDir, sprintf('%s_Biggest_Movers_Enrichment.png', clusterName))); close(f_movers);

    % =====================================================================
    % --- PLOT 3: METABOLIC DIET & STRATEGY ---
    if ~isempty(active_src_idx) || ~isempty(o2_idx)
        f_diet = figure('Name', 'Diet & Strategy', 'Position', [150, 150, 1500, 500], 'Visible', 'off'); tl_diet = tiledlayout(1, 2, 'TileSpacing', 'compact');
        
        ax1 = nexttile; hold on; grid on; diet_c = lines(length(active_src_idx)); 
        for s = 1:length(active_src_idx), plot(k_steps, src_traj_rel(s, :), 'LineWidth', 3, 'Color', diet_c(s,:), 'DisplayName', active_sources{s, 1}); end
        set(ax1, 'XDir', 'reverse'); xlabel(ax1, 'Surviving EFMs (K)'); ylabel(ax1, 'Abundance in EFM Pool (%)'); title(ax1, 'Carbon Substrate Abundance'); legend(ax1, 'Location', 'eastoutside');
        
        ax2 = nexttile; hold on; grid on;
        plot(k_steps, strat_aerobic, 'LineWidth', 3, 'Color', [0.00, 0.44, 0.74], 'DisplayName', 'Aerobic Respiration (+O2)');
        plot(k_steps, strat_ferment, 'LineWidth', 3, 'Color', [0.85, 0.32, 0.09], 'DisplayName', 'Fermentation/Anaerobic (-O2)');
        set(ax2, 'XDir', 'reverse'); xlabel(ax2, 'Surviving EFMs (K)'); ylabel(ax2, 'Abundance in EFM Pool (%)'); title(ax2, 'Energy Strategy: Respiration vs. Fermentation'); legend(ax2, 'Location', 'eastoutside');
        
        title(tl_diet, sprintf('Metabolic Strategy Shifts (%s)', upper(clusterName)), 'FontSize', 16, 'FontWeight', 'bold');
        saveas(f_diet, fullfile(plotDir, sprintf('%s_Metabolic_Strategy.png', clusterName))); close(f_diet);
    end

    % =====================================================================
    % --- PLOT: THE BIOLOGICAL SPECIALIZATION SIGNAL (FOLD-ENRICHMENT) ---
    % =====================================================================
    fprintf('  Extracting hidden specialization signals (Fold-Enrichment)...\n');
    
    % We need the relative % share of each subsystem at every K step
    sys_rel_share = zeros(length(subsys_names), length(k_steps));
    
    for t = 1:length(k_steps)
        step_surv = setdiff(cluster_initial_ids, removal_log.Removed_Global_ID(removal_log.Remaining_K > k_steps(t)));
        if isempty(step_surv), continue; end
        
        E_step_b = abs(E_full(:, step_surv)) > 1e-6;
        num_surv = length(step_surv);
        
        for s = 1:length(subsys_names)
            % What percentage of the SURVIVING EFMs use this subsystem?
            efms_using_sys = sum(any(E_step_b(subsys_map(subsys_names{s}), :), 1));
            sys_rel_share(s, t) = (efms_using_sys / num_surv) * 100;
        end
    end
    
    % Calculate Fold Enrichment: (Share at current K) / (Share at starting K)
    % Add a tiny epsilon to avoid division by zero for pathways that started dead
    fold_enrichment = zeros(length(subsys_names), length(k_steps));
    base_share = sys_rel_share(:, 1) + 1e-3; 
    
    for t = 1:length(k_steps)
        fold_enrichment(:, t) = sys_rel_share(:, t) ./ base_share;
    end
    
    % Find the truly "Special" pathways: those that MULTIPLIED the most
    % We filter out pathways that are just tiny noise (must reach at least 5% prevalence)
    final_share = sys_rel_share(:, end);
    valid_special = find(final_share > 5.0); 
    
    % Sort by maximum fold enrichment
    [~, sort_fold] = sort(fold_enrichment(valid_special, end), 'descend');
    
    % Grab the top 8 most highly enriched "special" pathways
    num_special = min(8, length(sort_fold));
    top_special_idx = valid_special(sort_fold(1:num_special));
    
    f_special = figure('Name', 'Biological Specialization', 'Position', [150, 150, 1200, 600], 'Visible', 'off');
    hold on; grid on;
    
    special_colors = turbo(num_special);
    
    for i = 1:num_special
        idx = top_special_idx(i);
        plot(k_steps, fold_enrichment(idx, :), 'LineWidth', 3, 'Color', special_colors(i,:), ...
            'DisplayName', sprintf('[%.1fx Spike] %s (Ends at %.1f%%)', fold_enrichment(idx, end), strrep(subsys_names{idx},'_',' '), final_share(idx)));
    end
    
    % Draw a baseline at 1.0 (No enrichment)
    yline(1.0, 'k--', 'LineWidth', 2, 'DisplayName', 'Baseline (No Enrichment)');
    
    set(gca, 'XDir', 'reverse'); 
    xlabel('Number of Surviving EFMs (K)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Fold-Enrichment (Multiplier vs. Initial Pool)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('The Specialization Signal: Identifying Unique Metabolic States (%s)', upper(clusterName)), 'FontSize', 14);
    legend('Location', 'eastoutside', 'FontSize', 10, 'Interpreter', 'none');
    
    saveas(f_special, fullfile(plotDir, sprintf('%s_Specialization_Fold_Enrichment.png', clusterName)));
    close(f_special);
end
fprintf('\n=======================================================\n');
fprintf('ANALYSIS COMPLETE: %s\n', datestr(now));
fprintf('LOG SAVED TO: %s\n', logFileName);
fprintf('=======================================================\n');
diary off; % STOP LOGGING