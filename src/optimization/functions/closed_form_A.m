function [A_opt, E, L, metrics] = closed_form_A(E, X, rxnNames, verbose, plotDir, clusterName, lambda_l21, cov_threshold, div_by_reactions, threshold_method)
    if nargin < 9, threshold_method = 'rowsum'; end
    fprintf('Method: %s\n', threshold_method);
    fprintf('=== Calculate A via closed form with REFIT ===\n');
    
    total_efms = size(E, 2);
    
    % --- 1. INITIAL SOLVE ---
    EtE = E' * E;
    [U, S, V] = svd(EtE, 'econ');
    s_vals = diag(S);
    s_inv = zeros(size(s_vals));
    keep_idx = s_vals > (1e-6 * max(s_vals));
    s_inv(keep_idx) = 1 ./ s_vals(keep_idx);
    EtE_pseudoinv = V * diag(s_inv) * U';
    
    A_initial = EtE_pseudoinv * E' * X * E * EtE_pseudoinv;
    A_initial = (A_initial + A_initial') * 0.5;
    
    % Baseline Diagnostic
    base_err = norm(X - (E * A_initial * E'), 'fro') / norm(X, 'fro');
    fprintf('  [Diagnostic] Baseline Error (Pre-Pruning): %.2f%%\n', base_err * 100);
    
    efm_lengths = sum(abs(E) > 1e-8, 1)';
    
    % --- 2. METRIC CALC & INITIAL PLOT ---
    if strcmpi(threshold_method, 'rowsum')
        fprintf("We are using rowsums \n")
        metric_values = sum(abs(A_initial - diag(diag(A_initial))), 2);
        metric_name = 'RowSums';
    else
        metric_values = diag(A_initial);
        metric_name = 'Diagonal';
    end
    if div_by_reactions
        metric_values = metric_values ./ efm_lengths;
    end 
    dynamic_threshold = max(metric_values) * cov_threshold; 
    if ~isempty(plotDir)
        plot_metric_distribution(metric_values, dynamic_threshold, metric_name, clusterName, plotDir, 'Pre_Refit_Distribution');
    end
    stage1_idx = find(metric_values > dynamic_threshold);
    stage1_dropped = find(metric_values <= dynamic_threshold);
    
    fprintf("  Stage 1 Pruning (L2): Keeping %d / %d EFMs\n", length(stage1_idx), total_efms);
    
    % --- INITIALIZE LOGGING ARRAYS (For Biological Script) ---
    log_Iter = [];
    log_K = [];
    log_GlobalID = [];
    
    current_K_tracker = total_efms;
    
    % Log Stage 1 Drops (Sort so worst are removed first at Iteration 0)
    if ~isempty(stage1_dropped)
        [~, sort_stage1] = sort(metric_values(stage1_dropped), 'ascend');
        stage1_dropped_sorted = stage1_dropped(sort_stage1);
        
        for i = 1:length(stage1_dropped_sorted)
            current_K_tracker = current_K_tracker - 1;
            log_Iter(end+1, 1) = 0; % Iteration 0 indicates pre-filter
            log_K(end+1, 1) = current_K_tracker;
            log_GlobalID(end+1, 1) = stage1_dropped_sorted(i);
        end
    end
    
    % Create the reduced E matrix for ADMM
    E_stage1 = E(:, stage1_idx);
    
    % =========================================================
    % STAGE 2: THE L2,1 SCALPEL (ADMM for exact row sparsity)
    % =========================================================
    
    % Pass the fraction, global indices, and current K directly to ADMM
    [~, admm_survivors_idx, admm_metrics, admm_log_Iter, admm_log_K, admm_log_GlobalID] = ...
        solve_A_L21_ADMM(E_stage1, X, 10, 2000, 1e-4, stage1_idx, current_K_tracker);
    
    % Append ADMM removal logs
    log_Iter = [log_Iter; admm_log_Iter];
    log_K = [log_K; admm_log_K];
    log_GlobalID = [log_GlobalID; admm_log_GlobalID];
    
    % Map the ADMM survivors back to the original full-size EFM indices
    final_selected_idx = stage1_idx(admm_survivors_idx);
    fprintf("  Stage 2 Pruning (ADMM): Keeping %d / %d EFMs\n", length(final_selected_idx), length(stage1_idx));
    
    if isempty(final_selected_idx)
        error('ADMM threshold too high! All EFMs were pruned.');
    end
    
    % --- WRITE THE REMOVAL LOG FOR BIOLOGICAL SCRIPT ---
    Optimal_K = length(final_selected_idx);
    log_OptimalK = repmat(Optimal_K, length(log_Iter), 1);
    log_Var = zeros(length(log_Iter), 1); % Biological script expects this column
    
    removal_table = table(log_Iter, log_K, log_GlobalID, log_Var, log_OptimalK, ...
        'VariableNames', {'Iteration', 'Remaining_K', 'Removed_Global_ID', 'Var_Explained', 'Optimal_K'});
    
    % Save to log_dir next to plotDir
    if ~isempty(plotDir)
        logDir = fullfile(fileparts(plotDir), 'log_dir');
        if ~exist(logDir, 'dir'), mkdir(logDir); end
        removal_log_path = fullfile(logDir, [clusterName '_removal_log.csv']);
        writetable(removal_table, removal_log_path);
        fprintf('  -> Saved complete removal log to: %s\n', removal_log_path);
    end

    % =========================================================
    % STAGE 3: THE UNPENALIZED REFIT (Recovering exact magnitudes)
    % =========================================================
    E_final = E(:, final_selected_idx); 
    EtE_final = E_final' * E_final;
    [Ur, Sr, Vr] = svd(EtE_final, 'econ');
    sr_vals = diag(Sr);
    sr_inv = zeros(size(sr_vals));
    keep_r = sr_vals > (1e-6 * max(sr_vals));
    sr_inv(keep_r) = 1 ./ sr_vals(keep_r);
    
    EtE_final_inv = Vr * diag(sr_inv) * Ur';
    A_refit = EtE_final_inv * E_final' * X * E_final * EtE_final_inv;
    A_refit = (A_refit + A_refit') * 0.5;
    
    % Map back to the full-sized optimal A matrix
    A_opt = zeros(size(A_initial));
    A_opt(final_selected_idx, final_selected_idx) = A_refit;
    
    % --- 4. POST-REFIT DIAGNOSTICS & PLOT ---
    X_recon_final = E * A_opt * E';
    final_err = norm(X - X_recon_final, 'fro') / norm(X, 'fro');
    fprintf('  [Diagnostic] Final Error after REFIT: %.2f%%\n', final_err * 100);
    
    % Optional: Plot the NEW distribution to see the energy consolidation
    if strcmpi(threshold_method, 'rowsum')
        new_metrics = sum(abs(A_opt - diag(diag(A_opt))), 2);
    else
        new_metrics = diag(A_opt);
    end
    if div_by_reactions
        new_metrics = new_metrics ./ efm_lengths;
    end 
    
    if ~isempty(plotDir)
        plot_metric_distribution(new_metrics, dynamic_threshold, metric_name, clusterName, plotDir, 'Post_Refit_Distribution');
    end
    
    % Final L matrix for sampling
    A_ridge = A_opt + eye(size(A_opt)) * 1e-8;
    L = chol(A_ridge, 'lower');
    
    reaction_variances = diag(X_recon_final);
    active_rxn_mask = reaction_variances > 1e-12;
    active_rxn_vars = reaction_variances(active_rxn_mask);
    
    if ~isempty(plotDir)
        f_rxn_var = figure('Visible', 'off');
        if ~isempty(active_rxn_vars)
            log_rxn_vars = log10(active_rxn_vars);
            histogram(log_rxn_vars, 50, 'Normalization', 'pdf', 'FaceColor', [0.8 0.4 0.2]);
            title(sprintf('Reconstructed Reaction Variances (%s)', clusterName), 'Interpreter', 'none');
            xlabel('log10(Reaction Variance)');
            ylabel('Probability Density');
            grid on;
            
            % Add basic stats to the plot
            stats_rxn = {
                sprintf('Total Rxns: %d', length(reaction_variances)), ...
                sprintf('Active Rxns (>1e-12): %d', length(active_rxn_vars)), ...
                sprintf('Median log10(Var): %.2f', median(log_rxn_vars))
            };
            annotation('textbox', [0.15, 0.75, 0.3, 0.15], 'String', stats_rxn, ...
                       'FitHeightToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'k');
        end
        saveas(f_rxn_var, fullfile(plotDir, sprintf('%s_Reaction_Variance_Density.png', clusterName)));
        close(f_rxn_var);
        
        rxn_participation_counts = sum(abs(E_final) > 1e-8, 2); 
        
        % Basic stats
        total_rxns = length(rxn_participation_counts);
        pruned_rxns = sum(rxn_participation_counts == 0);
        active_rxns = total_rxns - pruned_rxns;
        
        % Filter out zeros before log transformation
        active_counts = rxn_participation_counts(rxn_participation_counts > 0);
        
        f_participation = figure('Visible', 'off');
        
        if ~isempty(active_counts)
            % Transform the counts to log10 scale
            log_counts = log10(active_counts);
            
            % Plot the histogram
            histogram(log_counts, 50, 'FaceColor', [0.4 0.6 0.8]);
            
            title_str = sprintf('Log10 Reaction Participation in %d Selected EFMs (%s)', length(final_selected_idx), clusterName);
            title(title_str, 'Interpreter', 'none');
            
            % Explicit labels for the log scale
            xlabel('log10(Number of Participating EFMs)'); 
            ylabel('Number of Reactions'); 
            grid on;
            
            % Add a textbox to ensure we still see the zero-count (pruned) stats
            stats_str = {
                sprintf('Total Reactions: %d', total_rxns), ...
                sprintf('Pruned (0 EFMs): %d', pruned_rxns), ...
                sprintf('Active Reactions: %d', active_rxns)
            };
            annotation('textbox', [0.55, 0.75, 0.35, 0.15], 'String', stats_str, ...
                       'FitHeightToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'k');
        else
            title(sprintf('No Active Reactions Found (%s)', clusterName), 'Interpreter', 'none');
            text(0.5, 0.5, 'All reactions pruned', 'HorizontalAlignment', 'center');
        end
        % Save the plot
        saveas(f_participation, fullfile(plotDir, sprintf('%s_Log_Reaction_Participation_Hist.png', clusterName)));
        close(f_participation);
    else
        % Compute stats even if not plotting
        rxn_participation_counts = sum(abs(E_final) > 1e-8, 2); 
        total_rxns = length(rxn_participation_counts);
        pruned_rxns = sum(rxn_participation_counts == 0);
        active_rxns = total_rxns - pruned_rxns;
    end
    
    metrics = struct();
    % Update the metrics struct
    metrics.rxns_total = total_rxns;
    metrics.rxns_pruned = pruned_rxns;
    metrics.rxns_kept = active_rxns;
    
    metrics.base_error = base_err;
    metrics.final_error = final_err;
    metrics.efms_kept = length(final_selected_idx);
end