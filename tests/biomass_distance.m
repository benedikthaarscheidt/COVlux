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
usesecondmoment = config.params.second_moment_mat;
run_name_cluster = config.params.input_clustering_folder;
use_conditions   = isfield(config.params, 'use_conditions') && config.params.use_conditions;

%% 2. PATH DEFINITIONS
% COVlux Results (To get lost_rxns files)
if usebigbasis
    covluxBase = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_bigbasis');
else
    covluxBase = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_smallbasis');
end

clusteringBase = fullfile(projectRoot, config.paths.results_dir, 'Clustering');
if use_conditions
    covDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition');
    %expressionDir = '/Users/benedikthaarscheidt/M.Sc./master_thesis_second_moment/scripts/data_prep/clustered_data_v2/Ecoli_eBW4_nFeat2000_Dims7_Res0.6/data_files/grouped_by_condition';
else
    covDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files');
end

% Get all Run folders
d = dir(fullfile(covluxBase, 'Run_*')); 
d = d([d.isdir]);

if isempty(d)
    error('No Run folders found in %s', covluxBase);
end

% Extract folder names
folderNames = {d.name};

% Convert to datetime objects using the pattern in the folder names
runTimes = NaT(length(d), 1);
for i = 1:length(d)
    % Remove 'Run_' prefix
    timeStr = extractAfter(folderNames{i}, 'Run_');
    try
        % Parse: yyyy-MM-dd_HH-mm-ss
        runTimes(i) = datetime(timeStr, 'InputFormat', 'yyyy-MM-dd_HH-mm-ss');
    catch
        % Fallback to modification time if parsing fails
        runTimes(i) = datetime(d(i).datenum, 'ConvertFrom', 'datenum');
    end
end

% Sort by actual run time (newest first)
[~, sortIdx] = sort(runTimes, 'descend');

% Debug output

% Select LATEST by actual run time
resultsDir = fullfile(covluxBase, d(sortIdx(1)).name);


% Clustering Data (To get Covariance Matrices & Expression)
clusteringBase = fullfile(projectRoot, config.paths.results_dir, 'Clustering');
if use_conditions
    covInputDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition','mapped_for_benchmarking');
    %expressionDir = '/Users/benedikthaarscheidt/M.Sc./master_thesis_second_moment/scripts/data_prep/clustered_data_v2/Ecoli_eBW4_nFeat2000_Dims7_Res0.6/data_files/grouped_by_condition';
else
    covInputDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'MRAS_outputs');
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files','mapped_for_benchmarking');
end
if ~exist(covDir, 'dir'), covDir = fileparts(covDir); end % Fallback

% Models
modelPath = fullfile(projectRoot, config.paths.models_dir, config.model.model_file);

fprintf('Results Source:    %s\n', resultsDir);
fprintf('Covariance Source: %s\n', covDir);

%% 3. LOAD MODEL
fprintf('Loading Model...\n');
data = load(modelPath);
if isfield(data, 'pruned_ir'), model = data.pruned_ir; else, vars=fieldnames(data); model=data.(vars{1}); end
n_rxns = length(model.rxns);

% Identify Biomass
bio_idx = find(strcmp(model.rxns, targetBiomassName));
if isempty(bio_idx), bio_idx = find(contains(model.rxns, 'BIOMASS'),1); end
if isempty(bio_idx), error('Biomass reaction not found.'); end

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
    % Read lost reactions -> CORRECTED approach based on MASTER_COMPARISON_ULTIMATE
    try
        % Read file as lines instead of using import options
        raw_lines = readlines(fullfile(files(k).folder, filename));
        
        % Skip header if present (lines containing 'lost' or 'Var')
        if length(raw_lines) > 0 && (contains(raw_lines(1), 'lost') || contains(raw_lines(1), 'Var'))
            raw_lines(1) = [];
        end
        
        % Clean up and get unique lost reaction names
        covlux_names = strtrim(raw_lines(strlength(raw_lines) > 0));
        
        % Map to model indices
        [~, lost_covlux] = ismember(covlux_names, model.rxns);
        lost_idx = lost_covlux(lost_covlux > 0);
        
        % Get survivors (all reactions except lost ones)
        survivors_cov = setdiff(1:n_rxns, lost_idx)';
        lost_cov_count = length(lost_idx);
        
        % Create kept mask
        kept_cov = false(n_rxns, 1);
        kept_cov(survivors_cov) = true;
        
    catch ME
        fprintf('Error reading COVlux file for %s: %s\n', short_name, ME.message);
        kept_cov = true(n_rxns, 1); % Keep all if file empty/error
        lost_cov_count = 0;
    end
    
    % --- B. GET iMAT/FASTCORE SURVIVORS ---
    exprFile = fullfile(expressionDir,  [expr_short_name, '_mapped_to_model.csv']);
    RH = []; RL = [];
    
    if exist(exprFile, 'file')
        % 1. Load the pre-mapped file (Columns are ALREADY b-numbers)
        T = readtable(exprFile, 'VariableNamingRule', 'preserve');
        
        % 2. Extract Data directly (No intersect needed)
        % The MRAS script already handled operons and aliases
        gene_vals = mean(T{:, :}, 1, 'omitnan')';
        
        % 3. Create Structure for getExpressionSets
        valid_genes_struct = struct();
        valid_genes_struct.gene = T.Properties.VariableNames(:); % b-numbers
        valid_genes_struct.value = gene_vals(:);
        
        % 4. Run Selection
        % Now getting the correct ~4000 genes
        [RH, RL] = getExpressionSets(model, valid_genes_struct);
        
    else
        fprintf('Warning: Mapped expression file not found for %s\n', short_name);
        % Optional: Fallback to old behavior or skip
    end
    
    % iMAT
    idx_imat = run_imat_milp_standard(model, RH, RL);
    kept_imat = false(n_rxns, 1); kept_imat(idx_imat) = true;
    lost_imat_count = n_rxns - sum(kept_imat);
    
    % FASTCORE
    try
        fast_model = fastcore(model, RH, 1e-4);
        [is_kept, ~] = ismember(model.rxns, fast_model.rxns);
        kept_fast = is_kept;
        lost_fast_count = n_rxns - sum(kept_fast);
    catch
        kept_fast = false(n_rxns, 1);
        lost_fast_count = n_rxns;
    end
    
    % --- C. GAP ANALYSIS (Get Added Reactions) ---
    [gap_cov, added_cov]   = calculate_gap_and_patch(model, kept_cov, bio_idx);
    [gap_imat, added_imat] = calculate_gap_and_patch(model, kept_imat, bio_idx);
    [gap_fast, added_fast] = calculate_gap_and_patch(model, kept_fast, bio_idx);
    
    % --- D. VARIANCE ANALYSIS (On Fixed Models) ---
    % Load Covariance Matrix for this cluster
    
    if usesecondmoment
        covFile = fullfile(covDir, [cluster_full_name '_SecondMoment.csv']);
    else 
        covFile = fullfile(covDir, [cluster_full_name '_COV.csv']);
    end 
    % Try fallback name if full name fails
    if ~exist(covFile, 'file')
         covFile = fullfile(covDir, [short_name '_MRAS_COV.csv']);
    end
    
    var_cov = NaN; var_imat = NaN; var_fast = NaN;
    
    if exist(covFile, 'file')
        try
            T_cv = readtable(covFile, 'ReadRowNames', true, 'PreserveVariableNames', true);
            X = table2array(T_cv);
            rxns_X = string(T_cv.Properties.RowNames);
            
            % Map Model Rxns to X indices
            [~, map_Model_to_X] = ismember(model.rxns, rxns_X);
            
            % Calculate Total Variance in Data (Trace of X)
            total_var = trace(X);
            
            % Helper to calc captured variance
            calc_var = @(mask) sum(diag(X(map_Model_to_X(mask & map_Model_to_X>0), map_Model_to_X(mask & map_Model_to_X>0)))) / total_var * 100;
            
            % 1. COVlux Fixed
            fixed_cov = kept_cov; fixed_cov(added_cov) = true;
            var_cov = calc_var(fixed_cov);
            
            % 2. iMAT Fixed
            fixed_imat = kept_imat; fixed_imat(added_imat) = true;
            var_imat = calc_var(fixed_imat);
            
            % 3. Fastcore Fixed
            fixed_fast = kept_fast; fixed_fast(added_fast) = true;
            var_fast = calc_var(fixed_fast);
            
        catch
            % Covariance load failed
        end
    end
    
    % --- LOGGING ---
    fprintf('%-20s | %-5d %-5d %-5d | %-5.1f %-5.1f %-5.1f | %-5d %-5d %-5d\n', ...
        short_name, gap_cov, gap_fast, gap_imat, var_cov, var_fast, var_imat, ...
        lost_cov_count, lost_fast_count, lost_imat_count);
        
    Stats(k).Cluster = string(short_name);
    Stats(k).Gap_COV = gap_cov;
    Stats(k).Gap_FAST = gap_fast;
    Stats(k).Gap_iMAT = gap_imat;
    Stats(k).Var_COV = var_cov;
    Stats(k).Var_FAST = var_fast;
    Stats(k).Var_iMAT = var_imat;
    Stats(k).Lost_COV = lost_cov_count;
    Stats(k).Lost_FAST = lost_fast_count;
    Stats(k).Lost_iMAT = lost_imat_count;
end

%% 5. VISUALIZATION
outputPdfPath = fullfile(resultsDir, 'Functional_Gap_and_Variance.pdf');
T = struct2table(Stats);

fig = figure('Name', 'Gap vs Variance', 'Color', 'w', 'Position', [100 100 1200 1500]);
t = tiledlayout(3, 1, 'Padding', 'compact');

% Plot 1: Functional Gap
nexttile;
bar_data_gap = [T.Gap_COV, T.Gap_FAST, T.Gap_iMAT];
b1 = bar(bar_data_gap, 'grouped');
b1(1).FaceColor = [0 0.45 0.74]; b1(2).FaceColor = [0.93 0.69 0.13]; b1(3).FaceColor = [0.85 0.33 0.1];
ylabel('Gap Size (Reactions Added)');
title('1. Functional Gap (Distance to Biomass)');
legend({'COVlux', 'FASTCORE', 'iMAT'}, 'Location', 'best');
xticklabels(T.Cluster); xtickangle(45); grid on;

% Plot 2: Recovered Variance (Fixed Models)
nexttile;
bar_data_var = [T.Var_COV, T.Var_FAST, T.Var_iMAT];
b2 = bar(bar_data_var, 'grouped');
b2(1).FaceColor = [0 0.45 0.74]; b2(2).FaceColor = [0.93 0.69 0.13]; b2(3).FaceColor = [0.85 0.33 0.1];
ylabel('Explained Variance (%)');
title('2. Variance Captured by Fixed Functional Models');
xticklabels(T.Cluster); xtickangle(45); grid on;
ylim([0 100]);

% Plot 3: Lost Reactions (NEW)
nexttile;
bar_data_lost = [T.Lost_COV, T.Lost_FAST, T.Lost_iMAT];
b3 = bar(bar_data_lost, 'grouped');
b3(1).FaceColor = [0 0.45 0.74]; b3(2).FaceColor = [0.93 0.69 0.13]; b3(3).FaceColor = [0.85 0.33 0.1];
ylabel('Number of Lost Reactions');
title('3. Lost Reactions per Method');
legend({'COVlux', 'FASTCORE', 'iMAT'}, 'Location', 'best');
xticklabels(T.Cluster); xtickangle(45); grid on;

exportgraphics(fig, outputPdfPath);
fprintf('\nReport saved to: %s\n', outputPdfPath);

%% ========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function [gap_size, added_indices] = calculate_gap_and_patch(model, active_mask, bio_idx)
    % Finds min reactions to add to 'active_mask' to allow biomass > 0.1
    [m, n] = size(model.S);
    c = zeros(2*n, 1);
    
    % Penalty only for INACTIVE reactions
    inactive_idx = find(~active_mask);
    c(inactive_idx) = 1; c(inactive_idx + n) = 1; 
    
    Aeq = [model.S, -model.S]; beq = zeros(m, 1);
    A_ineq = zeros(1, 2*n); A_ineq(bio_idx) = -1; A_ineq(bio_idx+n) = 1;
    b_ineq = -0.1; % Target biomass flux
    
    lb_vec = zeros(2*n, 1); ub_vec = inf(2*n, 1);
    for i = 1:n
        if model.ub(i) > 0, ub_vec(i) = model.ub(i); else, ub_vec(i) = 0; end
        if model.lb(i) < 0, ub_vec(i+n) = abs(model.lb(i)); else, ub_vec(i+n) = 0; end
    end
    
    options = optimoptions('linprog', 'Display', 'none');
    try
        x = linprog(c, A_ineq, b_ineq, Aeq, beq, lb_vec, ub_vec, options);
        if isempty(x)
            gap_size = NaN; added_indices = [];
        else
            v_net = x(1:n) - x(n+1:end);
            % Which previously INACTIVE reactions are now carrying flux?
            flux_in_inactive = abs(v_net(inactive_idx));
            is_added = flux_in_inactive > 1e-5;
            
            gap_size = sum(is_added);
            added_indices = inactive_idx(is_added);
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