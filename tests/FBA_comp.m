% =========================================================================
% METHOD COMPARISON: MINIMAL MEDIUM FBA CAPACITY HEATMAPS
% =========================================================================
% 1. Parses COVLUX final lost reactions.
% 2. Runs iMAT and FASTCORE based on expression data.
% 3. Calculates BASE MODEL max production capacity on Minimal Media (100%).
% 4. Tests COVLUX, iMAT, and FASTCORE capacity on the same minimal media.
% 5. Generates a comparative Heatmap (Targets vs. Methods) per cluster.
% =========================================================================


%% 1. ROBUST PROJECT ROOT FINDER & INIT
if exist(fullfile(pwd, 'config', 'config.json'), 'file')
    projectRoot = pwd;
else
    currentSearchPath = fileparts(mfilename('fullpath')); projectRoot = '';
    while length(currentSearchPath) > 1 
        if exist(fullfile(currentSearchPath, 'config', 'config.json'), 'file')
            projectRoot = currentSearchPath; break;
        end
        newPath = fileparts(currentSearchPath);
        if strcmp(newPath, currentSearchPath), break; end; currentSearchPath = newPath;
    end
end
if isempty(projectRoot), error('Could not locate "config/config.json".'); end

configFile = fullfile(projectRoot, 'config', 'config.json');
config = jsondecode(fileread(configFile));

run_name_cluster = config.params.input_clustering_folder;
use_conditions   = isfield(config.params, 'use_conditions') && config.params.use_conditions;

if config.params.use_big_basis, covluxBase = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_bigbasis');
else, covluxBase = fullfile(projectRoot, config.paths.results_dir, 'COVlux_cov_smallbasis'); end

d = dir(fullfile(covluxBase, 'Run_*')); d = d([d.isdir]); [~, idx] = sort([d.datenum], 'descend');
resultsDir = fullfile(covluxBase, d(idx(1)).name);

clusteringBase = fullfile(projectRoot, config.paths.results_dir, 'Clustering');
if use_conditions
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files', 'grouped_by_condition','mapped_for_benchmarking');
else
    expressionDir = fullfile(clusteringBase, run_name_cluster, 'data_files','mapped_for_benchmarking');
end

%% 2. LOAD BASE MODEL & ESTABLISH MINIMAL MEDIUM
fprintf('Loading Base Model...\n');
modelPath = fullfile(projectRoot, config.paths.models_dir, config.model.model_file);
m_struct = load(modelPath);
if isfield(m_struct, 'pruned_ir'), model_base = m_struct.pruned_ir; else, fn=fieldnames(m_struct); model_base=m_struct.(fn{1}); end
[m_mets, n_rxns] = size(model_base.S);
rxnNames = model_base.rxns;

% Setup Base Bounds
if isfield(model_base, 'lb'), lb_base = model_base.lb; else, lb_base = zeros(n_rxns, 1); end
if isfield(model_base, 'ub'), ub_base = model_base.ub; else, ub_base = 1000 * ones(n_rxns, 1); end
if isfield(model_base, 'rev')
    for idx_r=1:n_rxns, if model_base.rev(idx_r)==0 && lb_base(idx_r)<0, lb_base(idx_r)=0; end, end
end

% --- ENFORCE TRUE BIOLOGICAL MINIMAL MEDIUM ---
uptake_idx = find(startsWith(lower(rxnNames), 'ex_') & endsWith(lower(rxnNames), '_b'));
lb_base(uptake_idx) = 0; ub_base(uptake_idx) = 0; 

minimal_substrates = {'ex_glc__d_e_b', 'ex_glc_e_b', 'ex_o2_e_b', 'ex_h2o_e_b', 'ex_h_e_b', ...
                      'ex_nh4_e_b', 'ex_pi_e_b', 'ex_so4_e_b', 'ex_k_e_b', 'ex_na1_e_b', ...
                      'ex_mg2_e_b', 'ex_ca2_e_b', 'ex_cl_e_b', 'ex_fe2_e_b', 'ex_fe3_e_b', ...
                      'ex_zn2_e_b', 'ex_mn2_e_b', 'ex_cu2_e_b', 'ex_cobalt2_e_b', ...
                      'ex_mobd_e_b', 'ex_ni2_e_b', 'ex_sel_e_b', 'ex_tungs_e_b', ...
                      'ex_thm_e_b', 'ex_cbl1_e_b', 'ex_niacin_e_b','ex_glyc_e_b', 'ex_pyr_e_b', 'ex_ade_e_b', 'ex_hxan_e_b', 'ex_ura_e_b','ex_leu__L_e_b'};
                      
for i = 1:length(minimal_substrates)
    idx_m = find(strcmpi(rxnNames, minimal_substrates{i}));
    if ~isempty(idx_m), ub_base(idx_m) = 1000; end
end
glc_idx = find(strcmpi(rxnNames, 'ex_glc__d_e_b') | strcmpi(rxnNames, 'ex_glc_e_b'));
if ~isempty(glc_idx), ub_base(glc_idx) = 10; end 
o2_idx = find(strcmpi(rxnNames, 'ex_o2_e_b'));
if ~isempty(o2_idx), ub_base(o2_idx) = 20; end

%% 3. FIND TARGETS & CALCULATE 100% BASELINE CAPACITY
target_rxns_mod = {};
for i = 1:length(rxnNames)
    r_name = lower(rxnNames{i});
    is_output = startsWith(r_name, 'ex_') || startsWith(r_name, 'dm_') || contains(r_name, 'biomass') || contains(r_name, 'atpm');
    has_kw = any(contains(r_name, {'trp','met','cys','glu','gln','ala','arg','asn','asp','his','ile','leu','lys','phe','pro','ser','thr','tyr','val','mta', 'atp', 'nadh', 'fad'}));
    if is_output && has_kw && ~endsWith(r_name, '_b'), target_rxns_mod{end+1} = rxnNames{i}; end
end
target_rxns_mod = unique(target_rxns_mod)';

fprintf('\nCalculating 100%% Baseline Capacity on Unpruned Model...\n');
options = optimoptions('linprog', 'Display', 'off'); 
base_capacity = zeros(length(target_rxns_mod), 1);

for p = 1:length(target_rxns_mod)
    t_idx = find(strcmp(rxnNames, target_rxns_mod{p}), 1);
    lb_test = lb_base; ub_test = ub_base;
    lb_test(t_idx) = 0; ub_test(t_idx) = 1000; % Protect exhaust pipe
    
    f = zeros(n_rxns, 1); f(t_idx) = -1;
    [~, fval, exitflag] = linprog(f, [], [], model_base.S, zeros(m_mets, 1), lb_test, ub_test, options);
    if exitflag == 1 && -fval > 1e-6, base_capacity(p) = -fval; end
end

% Filter out targets the base model can't even produce on minimal media
valid_mask = base_capacity > 0;
target_rxns_mod = target_rxns_mod(valid_mask);
base_capacity = base_capacity(valid_mask);
fprintf('Found %d valid targets capable of de novo synthesis in the base model.\n', length(target_rxns_mod));

%% 4. CLUSTER LOOP & METHOD COMPARISON
plotDir = fullfile(resultsDir, 'plots'); if ~exist(plotDir, 'dir'), mkdir(plotDir); end
files = dir(fullfile(resultsDir, '*_unique_lost_rxns_full.csv'));

method_names = {'COVLUX', 'iMAT', 'FASTCORE'};

for k = 1:length(files)
    filename = files(k).name;
    cluster_full_name = strrep(filename, '_unique_lost_rxns_full.csv', '');
    expr_name_parts = strsplit(cluster_full_name, '_MRAS');
    expr_short_name = expr_name_parts{1};
    
    fprintf('\n-------------------------------------------------------\n');
    fprintf('ANALYZING: %s\n', upper(expr_short_name));
    
    survivors = cell(1, 3);
    
    % --- A. COVLUX SURVIVORS ---
    try
        raw_lines = readlines(fullfile(files(k).folder, filename));
        if length(raw_lines)>0 && (contains(raw_lines(1), 'lost') || contains(raw_lines(1), 'Var')), raw_lines(1)=[]; end
        covlux_names = strtrim(raw_lines(strlength(raw_lines)>0));
        [~, lost_covlux] = ismember(covlux_names, model_base.rxns);
        lost_idx = lost_covlux(lost_covlux > 0);
        survivors{1} = setdiff(1:n_rxns, lost_idx)';
    catch, survivors{1} = []; end
    
    % --- B. EXPRESSION DATA -> iMAT & FASTCORE ---

    exprFile = fullfile(expressionDir, [expr_short_name, '_mapped_to_model.csv']);
    RH_indices = []; RL_indices = [];
    if exist(exprFile, 'file')
        T = readtable(exprFile, 'VariableNamingRule', 'preserve');
        file_genes = T.Properties.VariableNames;
        gene_vals  = mean(T{:, :}, 1, 'omitnan')';
        
        % FIX: Explicitly assign to a scalar struct to prevent struct array expansion
        valid_genes_struct = struct();
        valid_genes_struct.gene = file_genes(:); 
        valid_genes_struct.value = gene_vals(:);
        
        [RH_indices, RL_indices] = getExpressionSets(model_base, valid_genes_struct);
    end
    
    bio_idx = find(contains(model_base.rxns, 'BIOMASS'),1);
    if ~isempty(bio_idx), RH_indices = unique([RH_indices; bio_idx]); end
    
    fprintf('  -> Running iMAT...\n');
    try survivors{2} = run_imat_milp_strict(model_base, RH_indices, RL_indices); catch, survivors{2} = []; end
    
    fprintf('  -> Running FASTCORE...\n');
    try 
        kept_fastcore = fastcore(model_base, RH_indices, 1e-4); 
        [is_kept, ~] = ismember(model_base.rxns, kept_fastcore.rxns); 
        survivors{3} = find(is_kept); 
    catch, survivors{3} = []; end
    
    % --- C. TEST FBA CAPACITY PER METHOD ---
    method_results_perc = zeros(length(target_rxns_mod), 3);
    
    for m_idx = 1:3
        curr_surv = survivors{m_idx};
        if isempty(curr_surv), continue; end
        
        inactive_mask = true(n_rxns, 1); 
        inactive_mask(curr_surv) = false;
        
        for p = 1:length(target_rxns_mod)
            t_idx = find(strcmp(rxnNames, target_rxns_mod{p}), 1);
            
            lb_step = lb_base; ub_step = ub_base;
            lb_step(inactive_mask) = 0; ub_step(inactive_mask) = 0;
            
            % Protect the specific exhaust pipe
            lb_step(t_idx) = 0; ub_step(t_idx) = 1000; 
            
            f = zeros(n_rxns, 1); f(t_idx) = -1; 
            [~, fval, exitflag] = linprog(f, [], [], model_base.S, zeros(m_mets, 1), lb_step, ub_step, options);
            
            if exitflag == 1
                % Normalize against base capacity
                method_results_perc(p, m_idx) = (-fval / base_capacity(p)) * 100;
            end
        end
    end
    
    % --- D. PLOT COMPARATIVE HEATMAP ---
    f_comp = figure('Name', 'Method Comparison Capacity', 'Position', [200, 200, 600, max(500, 20*length(target_rxns_mod))], 'Visible', 'off');
    
    imagesc(method_results_perc); 
    colormap(parula); 
    cbar = colorbar; cbar.Label.String = 'Relative Max Production (%)';
    
    % X-Axis = Methods
    set(gca, 'XTick', 1:3, 'XTickLabel', method_names, 'FontSize', 10, 'FontWeight', 'bold');
    
    % Y-Axis = Targets
    set(gca, 'YTick', 1:length(target_rxns_mod), 'YTickLabel', target_rxns_mod, 'TickLabelInterpreter', 'none', 'FontSize', 8);
    title(sprintf('Minimal Media Capacity (%s)', upper(expr_short_name)));
    
    % Overlay values as text inside the heatmap for clarity
    for i = 1:length(target_rxns_mod)
        for j = 1:3
            val = method_results_perc(i, j);
            if val > 50, txt_col = 'k'; else, txt_col = 'w'; end
            text(j, i, sprintf('%.0f%%', val), 'HorizontalAlignment', 'center', 'Color', txt_col, 'FontSize', 8);
        end
    end
    
    saveas(f_comp, fullfile(plotDir, sprintf('%s_Method_Comparison_Heatmap.png', cluster_full_name)));
    close(f_comp);
    
    % Save CSV Data
    out_table = array2table(method_results_perc, 'VariableNames', method_names, 'RowNames', target_rxns_mod);
    writetable(out_table, fullfile(plotDir, sprintf('%s_Method_Comparison_Data.csv', cluster_full_name)), 'WriteRowNames', true);
end
fprintf('\nAll Method Comparison Heatmaps saved to: %s\n', plotDir);

%% HELPER FUNCTIONS (Preserved exactly as provided)
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
    if isempty(gpr_indices), RH_indices = []; RL_indices = []; return; end
    gpr_expr = rxnExpression(gpr_indices);
    high_thresh = prctile(gpr_expr, 75);
    low_thresh = prctile(gpr_expr, 25);
    RH_tmp = gpr_indices(gpr_expr >= high_thresh);
    RL_tmp = gpr_indices(gpr_expr <= low_thresh);
    RH_indices = RH_tmp(:);
    RL_indices = RL_tmp(:);
end