% =========================================================================
% EFM BASIS STRUCTURAL ANALYSIS (Raw Basis, No Redundancy Filtering)
% =========================================================================

%% 1. SETUP ENVIRONMENT & PATHS
currentScriptPath = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(currentScriptPath)); % src/optimization -> src -> root

configFile = fullfile(projectRoot, 'config', 'config.json');
if ~exist(configFile, 'file')
    error('Config file not found at: %s', configFile);
end
config = jsondecode(fileread(configFile));

usebigbasis      = config.params.use_big_basis;
SNR_THRESHOLD    = config.params.snr_threshold;
FLUX_NOISE_FLOOR = config.params.flux_noise_floor;

modelPath   = fullfile(projectRoot, config.paths.models_dir, config.model.model_file);
efmMatPath  = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{1});
efmMatPath2 = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{2});
efmMatPath3 = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{3});
efmMatPath4 = fullfile(projectRoot, config.paths.models_dir, config.model.efm_basis_files{4});

% Output Directory
outDir = fullfile(projectRoot, config.paths.results_dir, 'Raw_EFM_Analysis');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% 2. LOAD MODEL & EFM BASIS
fprintf('Loading Metabolic Model...\n');
data = load(modelPath);
if isfield(data, 'pruned_ir')
    model_ir = data.pruned_ir;
else
    vars = fieldnames(data);
    model_ir = data.(vars{1});
end
rxnNames = model_ir.rxns;

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

fprintf('Merging EFM Bases...\n');
if usebigbasis
    E_combined_raw = [[res1; flux1(idx1,:)], [res2; flux2(idx2,:)], [res3; flux3(idx3,:)], [res4; flux4(idx4,:)]];
else
    E_combined_raw = [[res1; flux1(idx1,:)], [res2; flux2(idx2,:)], [res4; flux4(idx4,:)]];
end

%% 3. FILTER NOISE / GHOSTS (Strictly NO Biomass or Jaccard Filtering)
fprintf('Running Strict SNR/Ghost Filter...\n');
[E_clean, ~] = quality_filter(E_combined_raw, SNR_THRESHOLD, FLUX_NOISE_FLOOR);
E_full = E_clean(2:end, :); % Strip residual row
n_efms = size(E_full, 2);
fprintf('Final Valid EFM Count for Analysis: %d\n', n_efms);

%% 4. ANALYSIS & PLOTTING: LENGTH & COVERAGE
E_bin = abs(E_full) > FLUX_NOISE_FLOOR;
efm_lengths = sum(E_bin, 1);
rxn_coverage = sum(E_bin, 2);

fig1 = figure('Name', 'EFM Length and Coverage', 'Position', [100, 100, 1200, 500], 'Color', 'w');

% Plot 1: EFM Lengths
subplot(1, 2, 1);
histogram(efm_lengths, 50, 'FaceColor', [0.2 0.6 0.5]);
title(sprintf('EFM Length Distribution (N=%d)', n_efms));
xlabel('Length (Number of Active Reactions)'); ylabel('Count'); grid on;
xline(median(efm_lengths), '--r', 'LineWidth', 1.5, 'Label', 'Median');

% Plot 2: Reaction Coverage
subplot(1, 2, 2);
histogram(rxn_coverage(rxn_coverage > 0), 50, 'FaceColor', [0.8 0.3 0.1]);
set(gca, 'YScale', 'log');
title('Reaction Coverage Frequency');
xlabel('Number of EFMs covering a reaction'); ylabel('Count (Log Scale)'); grid on;

saveas(fig1, fullfile(outDir, 'EFM_Distributions.png'));
close(fig1);

%% 5. ANALYSIS: SUBSYSTEM MAPPING
if isfield(model_ir, 'subSystems')
    subsys_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:length(rxnE)
        m_idx = find(strcmp(model_ir.rxns, rxnE(i)));
        if ~isempty(m_idx) && ~isempty(model_ir.subSystems{m_idx})
            sys = model_ir.subSystems{m_idx};
            if iscell(sys), sys = sys{1}; end
            if isempty(sys), sys = 'Unassigned'; end
            
            if ~isKey(subsys_map, sys)
                subsys_map(sys) = [];
            end
            subsys_map(sys) = [subsys_map(sys), i]; 
        end
    end
    
    sys_names = keys(subsys_map);
    sys_counts = zeros(length(sys_names), 1);
    
    for s = 1:length(sys_names)
        rxn_indices = subsys_map(sys_names{s});
        efms_touching = any(E_bin(rxn_indices, :), 1); 
        sys_counts(s) = sum(efms_touching);
    end
    
    [sorted_counts, s_idx] = sort(sys_counts, 'descend');
    top_n = min(50, length(sorted_counts)); 
    
    fig2 = figure('Name', 'Subsystem Coverage', 'Position', [150, 150, 1000, 600], 'Color', 'w');
    bar(sorted_counts(1:top_n), 'FaceColor', [0.3 0.5 0.7]);
    set(gca, 'XTick', 1:top_n, 'XTickLabel', sys_names(s_idx(1:top_n)), 'TickLabelInterpreter', 'none');
    xtickangle(45);
    title('Number of EFMs mapping to each Subsystem (Top 30)');
    ylabel('Number of EFMs'); grid on;
    
    saveas(fig2, fullfile(outDir, 'Subsystem_Mapping.png'));
    close(fig2);
end

%% 6. THE 50 SHORTEST EFMs TABLE
[sorted_lengths, sort_idx] = sort(efm_lengths, 'ascend');
num_to_extract = min(50, n_efms);

shortest_table_data = cell(num_to_extract, 3);
for i = 1:num_to_extract
    e_idx = sort_idx(i);
    active_rxns = rxnE(E_bin(:, e_idx));
    
    shortest_table_data{i, 1} = sprintf('EFM_%d', e_idx);
    shortest_table_data{i, 2} = sorted_lengths(i);
    shortest_table_data{i, 3} = strjoin(active_rxns, ', ');
end

T_short = cell2table(shortest_table_data, 'VariableNames', {'EFM_Original_Index', 'Length', 'Active_Reactions'});
writetable(T_short, fullfile(outDir, 'Top_50_Shortest_EFMs.csv'));
fprintf('Saved Top 50 Shortest EFMs to %s\n', fullfile(outDir, 'Top_50_Shortest_EFMs.csv'));


%% 5.5 IDENTIFY UNDERREPRESENTED SUBSYSTEMS & MISSING REACTIONS
fprintf('\n--- UNDERREPRESENTED SUBSYSTEMS & MISSING REACTIONS ---\n');

% 1. Find Subsystems with LOW coverage (e.g., < 20 EFMs)
low_cov_threshold = 50; 
low_sys_idx = find(sys_counts < low_cov_threshold);
[low_sys_counts, sort_low_idx] = sort(sys_counts(low_sys_idx), 'ascend');
low_sys_names = sys_names(low_sys_idx(sort_low_idx));

fprintf('Subsystems with fewer than %d EFMs (Consider re-enumerating these):\n', low_cov_threshold);
for s = 1:length(low_sys_names)
    fprintf('  - %-35s : %d EFMs\n', low_sys_names{s}, low_sys_counts(s));
end

% 2. Find specific Reactions with ZERO coverage
uncovered_idx = find(rxn_coverage == 0);
uncovered_rxns = rxnE(uncovered_idx);

if ~isempty(uncovered_idx)
    fprintf('\n Found %d reactions with ZERO EFM coverage. \n', length(uncovered_idx));
    
    uncovered_data = cell(length(uncovered_idx), 2);
    for i = 1:length(uncovered_idx)
        r_name = char(uncovered_rxns(i)); % Force char array
        
        % Map back to the original model to get the subsystem
        m_idx = find(strcmp(model_ir.rxns, r_name));
        if ~isempty(m_idx) && ~isempty(model_ir.subSystems{m_idx})
            sys = model_ir.subSystems{m_idx};
            if iscell(sys), sys = sys{1}; end
        else
            sys = 'Unassigned';
        end
        uncovered_data{i, 1} = r_name;
        uncovered_data{i, 2} = sys;
    end
    
    % Save to CSV grouped by Subsystem
    T_uncov = cell2table(uncovered_data, 'VariableNames', {'Reaction', 'SubSystem'});
    T_uncov = sortrows(T_uncov, 'SubSystem');
    uncov_file = fullfile(outDir, 'Zero_Coverage_Reactions.csv');
    writetable(T_uncov, uncov_file);
    fprintf('Saved explicit To-Do list of missing reactions to: %s\n', uncov_file);
else
    fprintf('\nGreat! All reactions in rxnE are covered by at least one EFM.\n');
end
fprintf('--------------------------------------------------------------\n');
%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [clean_mat, mask] = strict_ghost_filter(mat, snr, noise)
    % ONLY filters out EFMs violating the steady-state/SNR requirement.
    % NO redundancy/Jaccard filtering is applied here.
    
    residuals = mat(1, :);
    fluxes    = mat(2:end, :);
    n = size(mat, 2);
    mask = true(1, n);
    
    for k = 1:n
        res = residuals(k);
        v = abs(fluxes(:, k));
        
        active = v(v > noise);
        if isempty(active)
            min_flux = 0; 
        else
            min_flux = min(active); 
        end
        
        if res == 0
            ratio = Inf;
        elseif min_flux == 0
            ratio = 0;
        else
            ratio = min_flux / res;
        end
        
        if ratio < snr
            mask(k) = false;
        end
    end
    fprintf('Removed %d "Ghost" EFMs (Failed Steady-State/SNR).\n', sum(~mask));
    clean_mat = mat(:, mask);
end