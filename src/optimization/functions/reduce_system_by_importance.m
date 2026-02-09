function[E_out, X_out, mu_out, final_indices_relative, sort_idx_rxn] =
    reduce_system_by_importance(E_in, X_in, mu_in, use_weighted_qr,
                                apply_sorting) %
        REDUCE_SYSTEM_BY_IMPORTANCE %
        Pipeline : % 1. Sort Reactions by Activity(Mean MRAS).%
                   2. Aggressive Reduction by DATA RANK(rref on X).%
                   3. Cleanup Reduction by STOICHIOMETRIC RANK(rref on E).% %
                   Note : STRICTLY DOES NOT TOUCH EFM ORDER(COLUMNS).

                   %
                   -- -1. SORTING PHASE(ROWS ONLY)-- -
    if apply_sorting fprintf(
        '    [Step 1]: Sorting reactions by Activity (High->Low)...\n');
[ ~, sort_idx_rxn ] = sort(abs(mu_in), 'descend');

% Sort Rows(Reactions) E_sorted = E_in(sort_idx_rxn, :);
X_sorted = X_in(sort_idx_rxn, sort_idx_rxn);
mu_sorted = mu_in(sort_idx_rxn);
else fprintf(
    '    [Step 1]: Sorting Disabled. Using original reaction order.\n');
n_rxn = size(E_in, 1);
sort_idx_rxn = (1
                : n_rxn)';

    E_sorted = E_in;
X_sorted = X_in;
mu_sorted = mu_in;
end

        % -- -2. REDUCTION PHASE-- -
    if use_weighted_qr
        fprintf('    [Step 2]: Weighted QR on E (Aggressive Pruning)...\n');

variances = diag(X_sorted);
std_vals = sqrt(variances);
mean_abs = abs(mu_sorted);

        % --- FIX 1: Squared SNR for stronger separation ---
        snr_weights = (mean_abs ./ (std_vals + 1e-12));
        % this reduces the QR error on the reduced system like crazy since we
                only keep the informative reactions

                    snr_weights(snr_weights < 1e-8) = 1e-8;
        snr_weights(snr_weights > 1000) = 1000;

        W = diag(snr_weights);
        E_weighted_rows = W * E_sorted;

        [~, R_e, P_e] = qr(E_weighted_rows', 'vector');
        
        abs_diag_R = abs(diag(R_e));
        max_diag = max(abs_diag_R);
        
        % --- FIX 2: Aggressive Relative Tolerance (e.g. 0.1% of max signal) ---
        aggressive_tol = 1e-15 * max_diag; 
        
        rank_E = sum(abs_diag_R > aggressive_tol);
        
        final_indices_relative = P_e(1:rank_E);
        final_indices_relative = sort(final_indices_relative);
        
        fprintf('      -> Kept %d reactions (Aggressive Tol: %.2e, Max R: %.2e)\n', ...
            length(final_indices_relative), aggressive_tol, max_diag);
        
    else
        [~, indep_X] = rref(X_sorted);  
        E_temp = E_sorted(indep_X, :);
        [~, indep_E] = rref(E_temp');    
        
        final_indices_relative = indep_X(indep_E);
        
        fprintf('    -> Rank: %d (Data: %d -> Stoich: %d)\n', ...
            length(final_indices_relative), length(indep_X), length(indep_E));
    end

    % --- 3. BUILD OUTPUTS ---
    E_out = E_sorted(final_indices_relative, :);
    X_out = X_sorted(final_indices_relative, final_indices_relative);
    mu_out = mu_sorted(final_indices_relative);
end
