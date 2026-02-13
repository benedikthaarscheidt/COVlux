function [A_final, E_reduced, L_final, metrics] = precision_matrix(E, M, verbose, plotDir, clusterName)
% PRECISION-GUIDED SELECTION 
% this function needs improvement for safety and efficiency. Dont use in
% this state!
% 1. RANK: Prioritize "Unique" EFMs (High v'*Theta*v + High Expression).
% 2. PRUNE: Remove weakest EFMs one-by-one.
% 3. SELECT: Pick the model closest to the "Ideal Point" (0 Complexity, 100% Variance).
% 4. CALIBRATE: Use CVX to enforce non-negative weights on the final set.

    fprintf('\n=== Precision-Guided Selection (Closest to Ideal Point) ===\n');
    [n_rxn, n_efm] = size(E);

    % --- CONFIGURATION ---
    K_INIT          = 500;      % Start with Top 500
    LAMBDA_REG      = 1e-3;     % Precision stability
    BETA            = 0.5;      % 0.5 = Equal weight Independence vs Magnitude
    
    MIN_EFMS        = 5;        % Prune down to this number
    VAR_FLOOR       = 0;        % Absolute floor (safety)

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
    % 2. SCORES (Static Inputs)
    % =========================================================================
    % Expression Score (Magnitude)
    log_expr = log10(diag(M) + 1);
    log_expr = log_expr / max(log_expr);
    expr_raw = abs(E_norm)' * log_expr;
    expr_score = (expr_raw - min(expr_raw)) / (max(expr_raw) - min(expr_raw));
    
    % Precision Matrix (Independence)
    M_reg = M + LAMBDA_REG * eye(n_rxn);
    Theta = inv(M_reg);
    Theta = 0.5 * (Theta + Theta'); 

    % Precision Score
    prec_raw = diag(E_norm' * Theta * E_norm);
    prec_score = (prec_raw - min(prec_raw)) / (max(prec_raw) - min(prec_raw));

    % Combined Score & Sort
    final_score = (1 - BETA) * prec_score + BETA * expr_score;
    [~, sort_idx] = sort(final_score, 'descend');
    
    % Initial Pool
    initial_indices = sort_idx(1:min(K_INIT, n_valid));
    current_set = initial_indices;
    
    % =========================================================================
    % 3. BACKWARD ELIMINATION LOOP
    % =========================================================================
    fprintf('  Step 1: Pruning from %d down to %d...\n', length(current_set), MIN_EFMS);
    
    total_var_data = sum(diag(M));
    
    k_history = [];
    var_history = [];
    removed_history = [];
    
    iter = 0;
    while length(current_set) >= MIN_EFMS
        iter = iter + 1;
        k = length(current_set);
        E_cur = E_norm(:, current_set);
        
        % A. Precision Scores (Relative to current set)
        prec_sub = diag(E_cur' * Theta * E_cur);
        p_min = min(prec_sub); p_max = max(prec_sub);
        if p_max > p_min
            prec_norm = (prec_sub - p_min) / (p_max - p_min);
        else
            prec_norm = zeros(size(prec_sub));
        end
        
        % Combined Score
        expr_cur = expr_score(current_set);
        combined = (1 - BETA) * prec_norm + BETA * expr_cur;
        
        % B. Calculate Variance (
        try
            A_cur = solve_QR_regularized_new(E_cur, M, 1);
            A_cur = project_to_psd(A_cur);
            M_recon = E_cur * A_cur * E_cur';
            var_cur = sum(diag(M_recon)) / total_var_data;
        catch
            var_cur = 0;
        end
        
        % Store History
        k_history(end+1) = k;
        var_history(end+1) = var_cur;
        
        % D. Remove Weakest EFM
        [~, worst_local_idx] = min(combined);
        removed_efm = current_set(worst_local_idx);
        
        removed_history(end+1) = removed_efm;
        current_set(worst_local_idx) = []; 
        
        if mod(k, 50) == 0
            fprintf('    k=%3d | Var=%.1f%% | Removed %d\n', k, var_cur*100, valid_idx(removed_efm));
        end
    end
    
    % =========================================================================
    % 4. SELECTION: CLOSEST TO IDEAL POINT (0, 1)
    % =========================================================================
    fprintf('  Step 2: Selecting Model (Closest to Ideal Point)...\n');
    
    x = k_history(:);
    y = var_history(:);
    
    valid_mask = y >= VAR_FLOOR;
    
    if ~any(valid_mask)
        fprintf('    WARNING: No models met variance floor. Using Max Var.\n');
        [~, best_idx_hist] = max(y);
    else
        x_valid = x(valid_mask);
        y_valid = y(valid_mask);
        
        % Normalize Axes to [0, 1]
        % x (Complexity): 0 = Simple (few EFMs), 1 = Complex (many EFMs)
        x_norm = (x_valid - min(x_valid)) / (max(x_valid) - min(x_valid));
        
        % y (Variance): 0 = Poor, 1 = Perfect
        y_norm = (y_valid - min(y_valid)) / (max(y_valid) - min(y_valid));
        
        % Euclidean Distance to Ideal Point (x=0, y=1)
        dist_to_ideal = sqrt((x_norm - 0).^2 + (y_norm - 1).^2);
        
        % Find Minimum Distance
        [~, local_idx] = min(dist_to_ideal);
        
        valid_indices = find(valid_mask);
        best_idx_hist = valid_indices(local_idx);
    end
    
    best_k = x(best_idx_hist);
    best_var = y(best_idx_hist);
    
    % Reconstruct the Optimal Set
    start_size = length(initial_indices);
    n_removed = start_size - best_k;
    
    if n_removed > 0
        efms_to_remove = removed_history(1:n_removed);
        selected = setdiff(initial_indices, efms_to_remove, 'stable');
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
    
    % Use the Robust CVX Solver (Element-wise Non-Negative)
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
    
    fprintf('  --- Final Performance ---\n');
    fprintf('  Covariance Error: %.4f\n', cov_error);
    fprintf('  Variance Explained: %.2f%%\n', var_final*100);
    
    % Outputs
    final_selected_global = valid_idx(selected);
    E_reduced = E_sel; 
    A_final = zeros(n_efm, n_efm);
    
    % Map dense block to sparse matrix
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
        
        % Mark Selection
        plot(best_k, best_var, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        
        % Mark Ideal Direction (Top Left)
        if exist('valid_mask', 'var') && any(valid_mask)
             plot(min(x(valid_mask)), max(y(valid_mask)), 'g*', 'MarkerSize', 10);
        end
        
        title([clusterName ': Precision Selection (Closest to Ideal)']);
        xlabel('Number of EFMs (k)');
        ylabel('Variance Explained');
        set(gca, 'XDir', 'reverse'); 
        grid on;
        
        saveas(fig, fullfile(plotDir, [clusterName '_precision_selection.png']));
        close(fig);
    end
end