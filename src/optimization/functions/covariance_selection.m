function [A_final, E_reduced, L_final, metrics] = covariance_selection(E, M, verbose, plotDir, clusterName, BETA, SPARSITY_WEIGHT)
% COVARIANCE-GUIDED SELECTION (Geometric Batch Pruning)
%
% 1. RANK: Prioritize "High Signal" EFMs.
% 2. PRUNE: Remove weakest EFMs in batches (1% decay) for speed & stability.
% 3. SELECT: Pick model closest to Ideal Point (or use extracted metrics later).
% 4. CALIBRATE: Use CVX to enforce final weights.

    fprintf('\n=== Covariance Selection (Geometric Batch Pruning) ===\n');
    [n_rxn, n_efm] = size(E);
    
    % --- CONFIGURATION ---
    K_INIT          = 3000;      % Start with Top 3000
    MIN_EFMS        = 5;         % Prune down to this number
    VAR_FLOOR       = 0;         
    fprintf("Sparsity weight: %.2f \n", SPARSITY_WEIGHT);

    % =========================================================================
    % 1. PRE-PROCESSING
    % =========================================================================
    efm_norms = sqrt(sum(E.^2, 1));
    valid_efms = efm_norms > 1e-12;
    E_norm = zeros(size(E));
    E_norm(:, valid_efms) = E(:, valid_efms) ./ efm_norms(valid_efms);
    valid_idx = find(valid_efms);
    E_norm = E_norm(:, valid_idx);
    n_valid = length(valid_idx);

    % =========================================================================
    % 2. RANKING (Covariance + Magnitude)
    % =========================================================================
    % Covariance Score
    cov_all = diag(E_norm' * M * E_norm);
    
    % Expression Score
    log_expr = log10(diag(M) + 1);
    log_expr = log_expr / max(log_expr);
    expr_raw = abs(E_norm)' * log_expr;
    expr_score = (expr_raw - min(expr_raw)) / (max(expr_raw) - min(expr_raw));
    
    % Combined Score (Signal Strength)
    cov_score = (cov_all - min(cov_all)) / (max(cov_all) - min(cov_all));
    final_score = (1 - BETA) * cov_score + BETA * expr_score;
    
    [~, idx_sorted] = sort(final_score, 'descend');
    
    % Initial Pool
    initial_indices = idx_sorted(1:min(K_INIT, n_valid));
    current_set = initial_indices;
    
    % =========================================================================
    % 3. BACKWARD ELIMINATION LOOP (Batch Mode)
    % =========================================================================
    fprintf('  Step 1: Pruning from %d down to %d (1%% batch decay)...\n', length(current_set), MIN_EFMS);
    
    % Prepare Real-Time Log File
    removal_log_path = fullfile(plotDir, [clusterName '_removal_log.csv']);
    fid = fopen(removal_log_path, 'w');
    fprintf(fid, 'Iteration,Remaining_K,Removed_Global_ID,Var_Explained\n');
    fclose(fid);

    total_var_data = sum(diag(M));
    k_history = [];
    var_history = [];
    removed_history = [];
    
    iter_count = 0;
    
    while length(current_set) > MIN_EFMS
        iter_count = iter_count + 1;
        k = length(current_set);
        
        % Calculate Batch Size: 1% of current set, minimum 1
        n_to_remove = max(1, floor(k * 0.01)); 
        
        % Ensure we don't overshoot MIN_EFMS
        if (k - n_to_remove) < MIN_EFMS
            n_to_remove = k - MIN_EFMS;
        end
        
        E_cur = E_norm(:, current_set);
        
        % -------------------------------------------------------------
        % A. DYNAMIC SOLVE (Robust Ridge + PSD)
        % -------------------------------------------------------------
        try
            % 1. Stabilized Ridge Solve (Lambda=0.01)
            A_raw = solve_QR_regularized_fast(E_cur, M, 1e-2);
            
            % 2. Project to PSD 
            A_clean = project_to_psd(A_raw);
            
            % 3. Metrics
            current_variances = diag(A_clean);
            M_recon = E_cur * A_clean * E_cur';
            var_cur = sum(diag(M_recon)) / total_var_data;
            
            % Stats for logging
            min_val_A = min(A_clean(:));
            max_val_A = max(A_clean(:));
            
        catch
            warning('Solver failed at k=%d', k);
            current_variances = zeros(k, 1);
            var_cur = 0;
            min_val_A = 0; max_val_A = 0;
        end
        
        % -------------------------------------------------------------
        % B. SCORING (Max-Scaling to preserve floor)
        % -------------------------------------------------------------
        v_max = max(current_variances);
        if v_max > 1e-12
            cov_norm_dynamic = current_variances / v_max;
        else
            cov_norm_dynamic = zeros(k, 1);
        end
        
     
        expr_cur = expr_score(current_set);
        
        combined = (1 - BETA) * cov_norm_dynamic + BETA * expr_cur;
        
        % -------------------------------------------------------------
        % C. EXECUTION (Batch Cut)
        % -------------------------------------------------------------
        k_history(end+1) = k;
        var_history(end+1) = var_cur;
        
        % Identify Weakest EFMs in this batch
        [~, sort_idx] = sort(combined, 'ascend');
        worst_local_indices = sort_idx(1:n_to_remove);
        
        % Map to Global IDs
        batch_ids = current_set(worst_local_indices);
        
        % Save to history 
        removed_history = [removed_history; batch_ids(:)];
        
        % Write Batch to Log File
        fid = fopen(removal_log_path, 'a');
        for r_idx = batch_ids(:)'
            fprintf(fid, '%d,%d,%d,%.6f\n', iter_count, k, valid_idx(r_idx), var_cur);
        end
        fclose(fid);
        
        % Prune
        current_set(worst_local_indices) = []; 
        
        if mod(iter_count, 10) == 0 || n_to_remove == 1
            fprintf('    k=%4d | Var=%.1f%% | Batch=%d | Min A: %.2e | Max A: %.2e\n', ...
                k, var_cur*100, n_to_remove, min_val_A, max_val_A);
        end
    end
    
    % =========================================================================
    % 4. SELECTION: CLOSEST TO IDEAL POINT
    % =========================================================================
    fprintf('  Step 2: Selecting Model (Closest to Ideal Point)...\n');
    
    x = k_history(:);
    y = var_history(:);
    
    valid_mask = y >= VAR_FLOOR;
    
    if ~any(valid_mask)
        fprintf('    WARNING: No models met variance floor. Using Max Var.\n');
        [~, best_idx_hist] = max(y);
    else
        x_v = x(valid_mask);
        y_v = y(valid_mask);
        
        % Normalize [0, 1]
        x_n = (x_v - min(x_v)) / (max(x_v) - min(x_v) + 1e-12);
        y_n = (y_v - min(y_v)) / (max(y_v) - min(y_v) + 1e-12);
        
        % Distance to Ideal (0 Complexity, 100% Variance)
        dist_to_ideal = sqrt( SPARSITY_WEIGHT * (x_n - 0).^2 + (y_n - 1).^2 );
        [~, local_idx] = min(dist_to_ideal);
        
        valid_indices = find(valid_mask);
        best_idx_hist = valid_indices(local_idx);
    end
    
    best_k = x(best_idx_hist);
    best_var = y(best_idx_hist);
    
    % Reconstruct the Optimal Set using removed_history
    % We need to find how many EFMs were removed to reach best_k
    n_removed_total = length(initial_indices) - best_k;
    
    if n_removed_total > 0
        % Take the first N killed EFMs
        killed_efms = removed_history(1:n_removed_total);
        selected = setdiff(initial_indices, killed_efms, 'stable');
    else
        selected = initial_indices;
    end
    n_sel = length(selected);
    
    fprintf('\n  Optimal Cutoff: k=%d (Var=%.1f%%)\n', n_sel, best_var*100);
    
    % =========================================================================
    % 5. CALIBRATION (Rigorous CVX Optimization)
    % =========================================================================
    fprintf('  Step 3: Calibrating Final Weights (CVX)...\n');
    E_sel = E_norm(:, selected);
    
    A_full = solve_cvx(E_sel, M);
    
    % =========================================================================
    % 6. METRICS & OUTPUTS
    % =========================================================================
    M_recon_final = E_sel * A_full * E_sel';
    cov_error = norm(M - M_recon_final, 'fro') / norm(M, 'fro');
    var_final = sum(diag(M_recon_final)) / total_var_data;
    
    metrics.cov_error = cov_error;
    metrics.var_explained = var_final;
    metrics.n_selected = n_sel;
    metrics.k_history = k_history;
    metrics.var_history = var_history;
    
    % for Post-Hoc Analysis (Extracting other K values later)
    metrics.removed_history = removed_history;
    metrics.initial_indices = initial_indices;
    metrics.valid_idx = valid_idx;
    
    fprintf('  --- Final Performance ---\n');
    fprintf('  Covariance Error: %.4f\n', cov_error);
    fprintf('  Variance Explained: %.2f%%\n', var_final*100);
    
    % Outputs
    final_selected_global = valid_idx(selected);
    E_reduced = E_sel; 
    A_final = zeros(n_efm, n_efm);
    A_final(final_selected_global, final_selected_global) = A_full;
    
    try
        L_final = chol(A_full + 1e-12*eye(n_sel), 'lower');
    catch
        L_final = diag(sqrt(diag(A_full)));
    end
    
    % =========================================================================
    % 7. PLOT
    % =========================================================================
    if ~isempty(plotDir) && exist('plotDir', 'var')
        fig = figure('Visible', 'off');
        plot(x, y, '-b', 'LineWidth', 2); hold on;
        plot(best_k, best_var, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        
        if exist('valid_mask', 'var') && any(valid_mask)
             plot(min(x(valid_mask)), max(y(valid_mask)), 'g*', 'MarkerSize', 10);
        end
        
        title([clusterName ': Covariance Selection (Closest to Ideal)']);
        xlabel('Number of EFMs (k)');
        ylabel('Variance Explained');
        set(gca, 'XDir', 'reverse'); 
        grid on;
        
        saveas(fig, fullfile(plotDir, [clusterName '_covariance_selection.png']));
        close(fig);
    end
end