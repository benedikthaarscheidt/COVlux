%% MASTER_COMPARISON_ULTIMATE_ALL_PLOTS_EXTENDED
% REFACTORED: Uses config.json for paths.
% 1. Runs COVLUX, iMAT, FASTCORE.
% 2. WIRING: Covariance Analysis (Var, CorrM, Coh).
% 3. TOPOLOGY: Lost Components + DUAL CONNECTIVITY (LCC% & AlgConn).
% 4. INTEGRITY: AA% (Protein Synthesis) AND MTA% (Core Machinery).
% 5. CONSISTENCY: Falsepositives (Low Exp / High Flux).
% 6. REPORTING: Generates full CSV and Summary Table.
% =========================================================================

%% 1. CONFIGURATION
clear;

% -- -LOAD CONFIG-- - currentScriptPath = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(currentScriptPath));
% tests->root configFile = fullfile(projectRoot, 'config', 'config.json');
if
  ~exist(configFile, 'file'), error('Config file not found: %s', configFile);
end config = jsondecode(fileread(configFile));
resolve = @(p) fullfile(projectRoot, strrep(p, '/', filesep));

% -- -PATHS-- - modelPath =
    resolve(fullfile(config.paths.models_dir, config.model.model_file));

% Clustering Run Name(from config)
        clusteringRun = config.params.input_clustering_folder;
clusteringDir = fullfile(resolve(config.paths.results_dir), 'Clustering',
                         clusteringRun, 'data_files');
covInputDir = fullfile(clusteringDir, 'MRAS_outputs');
expressionDir = clusteringDir;

% COVLux Results Dir % Determine if big basis or
    small from config / folder check if config.params.use_big_basis covluxBase =
    fullfile(resolve(config.paths.results_dir), 'COVlux_cov_bigbasis');
else covluxBase =
    fullfile(resolve(config.paths.results_dir), 'COVlux_cov_smallbasis');
end

    % Find LATEST Run in COVLux results if not specified d =
    dir(fullfile(covluxBase, 'Run_*'));
d = d([d.isdir]);
if isempty (d)
  error('No COVLux run directories found in %s', covluxBase);
end[~, idx] = sort([d.datenum], 'descend');
runName = d(idx(1)).name;
resultsDir = fullfile(covluxBase, runName);

fprintf('Model: %s\n', modelPath);
fprintf('Clustering Data: %s\n', clusteringDir);
fprintf('COVLux Results: %s (Latest Run)\n', resultsDir);

targetBiomassName = 'BIOMASS_Ec_iML1515_WT_75p37M';

% -- -PHPP SETTINGS-- - n_points = 20;
glc_range = linspace(0.1, 20, n_points);
o2_range = linspace(0.1, 20, n_points);

% % 2. LOAD MODEL &PREPARE SUBSYSTEMS fprintf('Loading Base Model...\n');
m_struct = load(modelPath);
if isfield (m_struct, 'pruned_ir')
  , model = m_struct.pruned_ir;
else
  , fn = fieldnames(m_struct);
model = m_struct.(fn{1});
end[m, n] = size(model.S);
n_rxns_orig = n;
n_mets = m;

model_base = model;

% C.Run FBA model_test = changeObjective(model, targetBiomassName);
sol = optimizeCbModel(model_test, 'max');

if sol
  .f > 1e-6 fprintf('   [SUCCESS] The pruned model grows! Growth Rate: %.4f\n',
                    sol.f);
else
  fprintf('   [FAILURE] The pruned model DOES NOT grow (Rate: %.4e).\n', sol.f);
warning('Do not use this model. Check if essential reactions were pruned.');
end

        % -- -SUBSYSTEM CLEANUP-- -
    if isfield (model_base, 'subSystems') raw_subs = model_base.subSystems;
if ischar (raw_subs)
  , raw_subs = cellstr(raw_subs);
end flat_subs = strings(n, 1);
    for
      i = 1 : n if isempty (raw_subs{i}) flat_subs(i) = "Unassigned";
    elseif iscell(raw_subs{i}) flat_subs(i) = string(raw_subs{i} {1});
    else flat_subs(i) = string(raw_subs{i});
    end end else flat_subs = repmat("Unassigned", n, 1);
    end u_subs = unique(flat_subs);
    total_counts_map =
        containers.Map(cellstr(u_subs), zeros(length(u_subs), 1));
for
  i = 1 : length(u_subs) total_counts_map(char(u_subs(i))) =
              sum(flat_subs == u_subs(i));
end

    % Setup Media &Objective bio_idx =
    find(strcmp(model_base.rxns, targetBiomassName));
if isempty (bio_idx)
  , bio_idx = find(contains(model_base.rxns, 'BIOMASS'), 1);
end model_base.c = zeros(n, 1);
model_base.c(bio_idx) = 1;

% Setup PhPP Indices glc_idx = find(contains(model_base.rxns, 'EX_glc__D_e') &
                                        contains(model_base.rxns, '_b'),
                                    1);
o2_idx = find(
    contains(model_base.rxns, 'EX_o2_e') & contains(model_base.rxns, '_b'), 1);
do_phpp = ~isempty(glc_idx) && ~isempty(o2_idx);

% % 3. MAIN LOOP files =
    dir(fullfile(resultsDir, '*_unique_lost_rxns_full.csv'));
stats_out = struct([]);
opts_lp = optimoptions('linprog', 'Display', 'none');

fprintf(
    '\n%-28s | %-5s | %-5s | %-8s | %-6s | %-8s | %-6s | %-5s | %-6s | %-5s | %-5s | %-5s | %-5s | %-8s | %-8s | %-5s | %-5s | %-6s | %-30s\n',
    ... 'Cluster (Method)', 'Keep', 'Lost', 'Bio', 'LCC%', 'AlgConn', 'FVA%',
    'Block', 'falsepositives', 'PhPP', 'AA%', 'MTA%', 'DeadE', 'Var(X)',
    'CorrM', 'Coh.', 'nCC(L)', 'Sz(L)%', 'Most_Pruned_Cat');
fprintf('%s\n', repmat('-', 1, 240));

for
  k = 1 : length(files) model = model_base;
filename = files(k).name;
cluster_full_name = strrep(filename, '_unique_lost_rxns_full.csv', '');
expr_name_parts = strsplit(cluster_full_name, '_MRAS');
expr_short_name = expr_name_parts{1};

% -- -A.LOAD COVARIANCE-- - covFile =
    fullfile(covInputDir, [cluster_full_name '_COV.csv']);
if
  ~exist(covFile, 'file') % Fallback search d =
      dir(fullfile(covInputDir, [expr_short_name '*_COV*.csv']));
if
  ~isempty(d), covFile = fullfile(d(1).folder, d(1).name);
end end if exist (covFile, 'file') T_cov =
    readtable(covFile, 'ReadRowNames', true, 'PreserveVariableNames', true);
rxnNames_X = string(T_cov.Properties.RowNames);
X_full = table2array(T_cov);
[ ~, map_X_to_Model ] = ismember(rxnNames_X, model.rxns);
has_covariance = true;
else has_covariance = false;
end

            % -- -B.PREPARE EXPRESSION-- -
        % Assume format is cluster_X.csv or
    similar in dataDir % Try standard naming from clustering output exprFile =
    fullfile(expressionDir, [ expr_short_name, '.csv' ]);

RH_indices = [];
RL_indices = [];
current_RL_set = [];

% Note : If expression file not found exactly,
    try logCPM pattern if ~exist(exprFile, 'file') exprFile =
        fullfile(expressionDir, [ expr_short_name, '_logCPM.csv' ]);
end

    if exist (exprFile, 'file') T =
        readtable(exprFile, 'VariableNamingRule', 'preserve');
        % Filter for valid gene columns (usually columns after metadata)
        % Heuristic: intersect with model genes
        [valid_genes, ~, ~] = intersect(T.Properties.VariableNames, model.genes);

        if
          ~isempty(valid_genes) gene_vals = mean(
              T{ :, valid_genes}, 1, 'omitnan')'; [RH_indices, RL_indices] =
              getExpressionSets(model, struct('gene', {valid_genes}, 'value',
                                              gene_vals));
        current_RL_set = RL_indices;
        end end if ~isempty(bio_idx) RH_indices = unique([RH_indices; bio_idx]);
        end

                % -- -C.GET SURVIVORS-- -
            survivors = cell(1, 3);
        try % COVLUX Survivors(Inverse of Lost) raw_lines =
            readlines(fullfile(files(k).folder, filename));
        if length (raw_lines)
          > 0 &&
              (contains(raw_lines(1), 'lost') || contains(raw_lines(1), 'Var')),
              raw_lines(1) = [];
        end covlux_names = strtrim(raw_lines(strlength(raw_lines) > 0));
        [ ~, lost_covlux ] = ismember(covlux_names, model.rxns);
        lost_idx = lost_covlux(lost_covlux > 0);
        survivors{1} = setdiff(1 : n, lost_idx)'; catch, survivors{1} = [];
        end

            % iMAT and FASTCORE try survivors{2} =
            run_imat_milp_strict(model, RH_indices, RL_indices);
        catch, survivors{2} = [];
        end try kept_fastcore = fastcore(model, RH_indices, 1e-4);
        [ is_kept, ~] = ismember(model.rxns, kept_fastcore.rxns);
        survivors{3} = find(is_kept);
        catch, survivors{3} = [];
        end

                % -- -D.RUN ANALYSIS-- -
            names = {'COVLUX', 'iMAT', 'FASTCORE'};
    for
      m_idx = 1 : 3 model = model_base;
    curr_survivors_idx = survivors{m_idx};

    if
      ~ismember(bio_idx, curr_survivors_idx) %
          fprintf('   WARNING: Biomass reaction is in LOST set for %s!\n',
                  names{m_idx});
    curr_survivors_idx = union(curr_survivors_idx, bio_idx);
    lost_indices = setdiff(1 : n_rxns_orig, curr_survivors_idx);
    end curr_name = names{m_idx};
    num_kept = length(curr_survivors_idx);
    num_lost = n_rxns_orig - num_kept;
    lost_indices = setdiff(1 : n_rxns_orig, curr_survivors_idx);

    % ...(METRICS CALCULATION - Same logic as before)... %
        TOPOLOGY OF LOST num_lost_components = 0;
    largest_lost_fraction = 0;
    if num_lost
      > 0 S_lost = model.S( :, lost_indices);
            Adj_lost = [sparse(n_mets, n_mets), abs(S_lost) ~= 0; (abs(S_lost) ~= 0)', sparse(num_lost, num_lost)];
            [bins, ~] = conncomp(graph(Adj_lost));
            rxn_bins = bins(n_mets+1 : end);
            unique_bins = unique(rxn_bins);
            num_lost_components = length(unique_bins);
            if num_lost_components > 0
               [~,~,bin_idx] = unique(rxn_bins); bin_counts = accumarray(bin_idx, 1);
               largest_lost_fraction = (max(bin_counts) / num_lost) * 100;
            end
        end
        
        % SUBSYSTEM
        most_pruned_sub = "None"; most_pruned_pct = 0;
        if ~isempty(lost_indices)
            lost_subs = flat_subs(lost_indices); u_lost_subs = unique(lost_subs); max_pct = -1;
            for s = 1:length(u_lost_subs)
                this_sub = u_lost_subs(s);
                if this_sub == "Unassigned" || this_sub == "" || contains(this_sub, "Exchange", 'IgnoreCase', true), continue; end
                count_lost = sum(lost_subs == this_sub); count_total = total_counts_map(char(this_sub));
                if count_total > 5
                    pct_lost = count_lost / count_total;
                    if pct_lost > max_pct, max_pct = pct_lost; most_pruned_sub = this_sub; end
                end
            end
            if max_pct > -1, most_pruned_pct = max_pct; end
        end
        
        % WIRING
        var_per_rxn=0; coherence=0; corr_mass=0;
        if has_covariance && num_kept > 1
            keep_X_mask = ismember(map_X_to_Model, curr_survivors_idx);
            if sum(keep_X_mask) > 1 
                X_sub = X_full(keep_X_mask, keep_X_mask);
                var_total = trace(X_sub); var_per_rxn = var_total / size(X_sub, 1);
                off_diag = X_sub - diag(diag(X_sub)); corr_mass = sum(abs(off_diag(:)));
                coherence = corr_mass / var_total;
            end
        end
        
        % METABOLIC & GRAPH TOPOLOGY
        model_test = model;
        if ~isempty(lost_indices), model_test.lb(lost_indices)=0; model_test.ub(lost_indices)=0; end
        
        [~, fval, flag] = linprog(-model_test.c, [], [], model_test.S, zeros(n_mets,1), model_test.lb, model_test.ub, opts_lp);
        max_biomass = (flag == 1) * -fval;
        
        blocked_count = 0; falsepositives_count = 0;
        if num_kept > 0
            tol_fva = 1e-6;
            for i = 1:num_kept
                rxn_id = curr_survivors_idx(i); c_fva = zeros(n_rxns_orig, 1); c_fva(rxn_id) = -1;
                [~, f_max, flag_fva] = linprog(c_fva, [], [], model_test.S, zeros(n_mets,1), model_test.lb, model_test.ub, opts_lp);
                v_max_val = abs((flag_fva==1) * f_max);
                if v_max_val < tol_fva, blocked_count = blocked_count + 1;
                elseif ismember(rxn_id, current_RL_set), falsepositives_count = falsepositives_count + 1; end
            end
            functional_ratio = (num_kept - blocked_count) / num_kept;
        else, functional_ratio = 0; end
        
        S_pruned = model.S; S_pruned(:, lost_indices) = 0; 
        Adj = [sparse(n_mets, n_mets), abs(S_pruned) ~= 0; (abs(S_pruned) ~= 0)', sparse(n_rxns_orig, n_rxns_orig)];
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
        
        avg_phpp = 0;
        if do_phpp && max_biomass > 1e-6
             phpp_vals = [];
             for i_g = 1:n_points
                 for j_o = 1:n_points
                     model_p = model_test;
                     model_p.ub(glc_idx) = glc_range(i_g); model_p.ub(o2_idx) = o2_range(j_o);
                     [~, f_p, fl_p] = linprog(-model_p.c, [], [], model_p.S, zeros(n_mets,1), model_p.lb, model_p.ub, opts_lp);
                     if fl_p == 1, phpp_vals(end+1) = -f_p; end
                 end
             end
             if ~isempty(phpp_vals), avg_phpp = mean(phpp_vals(phpp_vals > 0.01)); end
        end
        
        % INTEGRITY
        aa_tasks = {'ala__L_c','arg__L_c','asn__L_c','asp__L_c','cys__L_c','gln__L_c','glu__L_c','gly_c','his__L_c','ile__L_c','leu__L_c','lys__L_c','met__L_c','phe__L_c','pro__L_c','ser__L_c','thr__L_c','trp__L_c','tyr__L_c','val__L_c'};
        passed_aa = 0;
        for ti = 1:length(aa_tasks)
            metIdx = find(strcmp(model.mets, aa_tasks{ti}));
            if ~isempty(metIdx)
                model_check = model_test; 
                model_check.S = [model_check.S, sparse(n_mets, 1)]; model_check.S(metIdx, end) = -1;
                model_check.c = [zeros(n,1); 1]; model_check.lb = [model_check.lb; 0]; model_check.ub = [model_check.ub; 1000];
                [~, f_aa, fl_aa] = linprog(-model_check.c, [], [], model_check.S, zeros(n_mets,1), model_check.lb, model_check.ub, opts_lp);
                if fl_aa == 1 && -f_aa > 1e-3, passed_aa = passed_aa + 1; end
            end
        end
        aa_pct = (passed_aa / length(aa_tasks)) * 100;
        
        mta_tasks = {'atp_c', 'ctp_c', 'gtp_c', 'utp_c', 'datp_c', 'dctp_c', 'dgtp_c', 'dttp_c', 'nad_c', 'nadh_c', 'nadp_c', 'nadph_c', 'fad_c', 'coa_c', 'pheme_c', 'pydx5p_c', 'pyr_c', 'accoa_c', 'succoa_c', 'akg_c', 'oaa_c'};
        passed_mta = 0;
        for ti = 1:length(mta_tasks)
            metIdx = find(strcmp(model.mets, mta_tasks{ti}));
            if ~isempty(metIdx)
                model_check = model_test; 
                model_check.S = [model_check.S, sparse(n_mets, 1)]; model_check.S(metIdx, end) = -1;
                model_check.c = [zeros(n,1); 1]; model_check.lb = [model_check.lb; 0]; model_check.ub = [model_check.ub; 1000];
                [~, f_mta, fl_mta] = linprog(-model_check.c, [], [], model_check.S, zeros(n_mets,1), model_check.lb, model_check.ub, opts_lp);
                if fl_mta == 1 && -f_mta > 1e-3, passed_mta = passed_mta + 1; end
            end
        end
        mta_pct = (passed_mta / length(mta_tasks)) * 100;
        
        has_prod = any(S_pruned > 0, 2); has_cons = any(S_pruned < 0, 2);
        internal_mets = ~contains(model.mets, '_e') & ~contains(model.mets, '_b');
        dead_ends = internal_mets & ((has_prod & ~has_cons) | (~has_prod & has_cons));
        num_dead_ends = sum(dead_ends);
        
        % PRINT & SAVE
        sub_str = sprintf('%s (%.3f)', most_pruned_sub, most_pruned_pct);
        if strlength(sub_str) > 30, sub_str = extractBefore(sub_str, 30); end
        row_name = sprintf('%s (%s)', expr_short_name, curr_name);
        if strlength(row_name) > 28, row_name = extractBefore(row_name, 28); end
        
        fprintf('%-28s | %-5d | %-5d | %-8.4f | %-6.1f | %-8.4f | %-6.1f | %-5d | %-6d | %-5.4f | %-5.0f | %-5.0f | %-5d | %-8.2e | %-8.2e | %-5.2f | %-5d | %-6.1f | %-30s\n', ...
            row_name, num_kept, num_lost, max_biomass, lcc_pct, alg_conn, functional_ratio*100, ...
            blocked_count, falsepositives_count, avg_phpp, aa_pct, mta_pct, num_dead_ends, var_per_rxn, corr_mass, coherence, ...
            num_lost_components, largest_lost_fraction, sub_str);
            
        s_idx = length(stats_out)+1;
        stats_out(s_idx).Cluster = expr_short_name;
        stats_out(s_idx).Method = curr_name;
        stats_out(s_idx).KeptCount = num_kept;
        stats_out(s_idx).LostCount = num_lost;
        stats_out(s_idx).Biomass = max_biomass;
        stats_out(s_idx).LCC_Pct = lcc_pct;                
        stats_out(s_idx).Alg_Connectivity = alg_conn;      
        stats_out(s_idx).FunctionalPct_FVA = functional_ratio * 100;
        stats_out(s_idx).BlockedSurvivors = blocked_count;
        stats_out(s_idx).Falsepositives = falsepositives_count;
        stats_out(s_idx).PhPP_Avg = avg_phpp;
        stats_out(s_idx).AA_Synthesis_Pct = aa_pct;
        stats_out(s_idx).MTA_Core_Pct = mta_pct;
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
writetable(struct2table(stats_out), fullfile(resultsDir, 'FULL_METHOD_COMPARISON_WITH_DUAL_CONN.csv'));
fprintf('\nDone. Saved to FULL_METHOD_COMPARISON_WITH_DUAL_CONN.csv\n');

%% 7. PROFESSOR SUMMARY REPORT
if ~isempty(stats_out)
    T_results = struct2table(stats_out);
    T_results.Failed = T_results.Biomass < 1e-6; 
    summary = grpstats(T_results, 'Method', {'mean', 'std', 'sum'}, ...
        'DataVars', {'Biomass', 'AA_Synthesis_Pct', 'MTA_Core_Pct', 'Falsepositives', 'LCC_Pct', 'Alg_Connectivity', 'Failed'});
    prof_report = table();
    prof_report.Method = summary.Method;
    prof_report.Success_Rate = 100 * (1 - (summary.sum_Failed ./ summary.GroupCount));
    prof_report.Avg_Integrity_AA = summary.mean_AA_Synthesis_Pct;
    prof_report.Avg_Core_MTA = summary.mean_MTA_Core_Pct;
    prof_report.Avg_Connectivity_LCC = summary.mean_LCC_Pct;
    prof_report.Avg_Robustness_AlgConn = summary.mean_Alg_Connectivity;
    prof_report.Avg_Falsepositives = summary.mean_Falsepositives;
    fprintf('\n\n=== SUMMARY REPORT FOR PROFESSOR ===\n');
    disp(prof_report);
    writetable(prof_report, fullfile(resultsDir, 'PROFESSOR_SUMMARY_REPORT.csv'));
end

%% FUNCTION: iMAT MILP
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
    % Uses mapExpressionToReactions (requires COBRA Toolbox)
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