function [A_final, E_reduced, L_final, metrics] = solve_efm_hybrid_selection_2(E, M, verbose, plotDir, clusterName)
% HYBRID SELECTION – Covariance target, variance‑driven forward selection
%   Then compute optimal full covariance A via regularized pseudoinverse.
%
%   - Normalize EFMs → flux ratios only, no magnitude bias.
%   - Score = v' * M * v (covariance explained) – orders EFMs by individual utility.
%   - Forward selection: add ALL EFMs that improve variance explained.
%   - Stop only at safety cap (MAX_EFMS) or when candidates exhausted.
%   - Compute A = pinv(E_norm) * M * pinv(E_norm)' with truncated SVD.
%   - Output: normalized EFMs + full covariance matrix A.


    fprintf('\n=== Hybrid EFM Selection (Full Covariance A) ===\n');
    [n_rxn, n_efm] = size(E);

    
    MAX_EFMS            = 500;    
    SVD_TOL             = 1e-6;    
    % ---------------------------------------------------------------------
    % 1. NORMALIZE EFMs 
    % ---------------------------------------------------------------------
    efm_norms = sqrt(sum(E.^2, 1));
    valid_efms = efm_norms > 1e-12;
    E_norm = zeros(size(E));
    E_norm(:, valid_efms) = E(:, valid_efms) ./ efm_norms(valid_efms);
    valid_idx = find(valid_efms);
    E_norm = E_norm(:, valid_idx);
    n_valid = length(valid_idx);

    % ---------------------------------------------------------------------
    % 2. SCORE EFMs – v' * M * v 
    % ---------------------------------------------------------------------
    fprintf('  Scoring %d EFMs by v'' M v...\n', n_valid);
    scores = diag(E_norm' * M * E_norm);
    [scores_sorted, idx_sorted] = sort(scores, 'descend');

    % ---------------------------------------------------------------------
    % 3. FORWARD SELECTION – ADD ALL EFMs THAT IMPROVE VARIANCE EXPLAINED
    % ---------------------------------------------------------------------
    fprintf('\n  Step 1: Forward selection (combined variance+covariance score)...\n');

    selected = [];
    M_recon = zeros(n_rxn);
    total_var = sum(diag(M));
    var_explained = 0;
    
    alpha = 0.5;   
    
    for k = 1:length(idx_sorted)
        if length(selected) >= MAX_EFMS
            fprintf('    Stopping: reached max EFMs (%d).\n', MAX_EFMS);
            break;
        end
    
        efm_local_idx = idx_sorted(k);
        v = E_norm(:, efm_local_idx);
    
        residual = M - M_recon;
    
        % --- Variance improvement ---
        residual_var = diag(residual);
        var_improve = max(0, sum(v.^2 .* residual_var));
    
        % --- Covariance improvement ---
        residual_cov = residual;
        residual_cov(1:n_rxn+1:end) = 0;   % zero diagonal
        cov_improve = max(0, v' * residual_cov * v);
    
        % --- Combined improvement ---
        combined_improve = (1-alpha) * var_improve + alpha * cov_improve;
    
        if combined_improve > 1e-12
            weight = max(0, v' * residual * v);   % full M
            M_recon = M_recon + weight * (v * v');
    
            new_var = sum(diag(M_recon)) / total_var;
            gain = new_var - var_explained;
    
            if gain > 1e-15   % keep only if actually improves variance (numerical safety)
                selected = [selected, efm_local_idx];
                var_explained = new_var;
    
                if mod(k, 10) == 0 || k == 1
                    fprintf('    k=%4d: added EFM %5d | var = %.2f%% | gain = %.4f%% | combined = %.2e\n', ...
                        k, valid_idx(efm_local_idx), var_explained*100, gain*100, combined_improve);
                end
            end
        end
    end
    fprintf('\n  Forward selection complete: %d EFMs selected, %.1f%% variance explained.\n', ...
        length(selected), var_explained*100);

     % ---------------------------------------------------------------------
    % 4. COMPUTE OPTIMAL A via TRUNCATED SVD (REGULARIZED PSEUDOINVERSE)
    %    ***** CRITICAL: USE NORMALIZED EFMs *****
    % ---------------------------------------------------------------------
    fprintf('\n  Step 2: Computing optimal A via truncated SVD...\n');
    E_sel = E_norm(:, selected);          
    n_sel = size(E_sel, 2);

    A_full = solve_QR_regularized_new(E_sel, M, 1);
    A_full = project_to_psd(A_full);
    % ---------------------------------------------------------------------
    % 5. METRICS 
    % ---------------------------------------------------------------------
    M_recon_final = E_sel * A_full * E_sel';
    cov_error = norm(M - M_recon_final, 'fro') / norm(M, 'fro');
    var_final = sum(diag(M_recon_final)) / total_var;

    % --- Correlation pattern match (off‑diagonal) ---
    dvec = sqrt(diag(M)); dvec(dvec<1e-12)=1;
    D_inv = diag(1./dvec);
    C_data = D_inv * M * D_inv; C_data(1:n_rxn+1:end)=0;
    d_recon = sqrt(diag(M_recon_final)); d_recon(d_recon<1e-12)=1;
    D_recon_inv = diag(1./d_recon);
    C_recon = D_recon_inv * M_recon_final * D_recon_inv; C_recon(1:n_rxn+1:end)=0;
    pattern_correlation = corr(C_data(:), C_recon(:));

    
    metrics = struct();
    metrics.cov_error = cov_error;
    metrics.var_explained = var_final;
    metrics.pattern_correlation = pattern_correlation;
    metrics.n_selected = length(selected);
    metrics.n_active = n_sel;
  
    fprintf('\n  --- Final Performance (Stable Truncated SVD) ---\n');
    fprintf('  Covariance reconstruction error: %.4f\n', cov_error);
    fprintf('  Variance explained: %.2f%%\n', var_final*100);
    fprintf('  Correlation pattern match (off‑diag): %.4f\n', pattern_correlation);
    fprintf('  Active EFMs: %d (all selected are kept)\n', n_sel);
    fprintf('  Min value in A: %.6e\n', min(A_full(:)));

    
    final_selected_global = valid_idx(selected);
    E_reduced = E_sel;                     
    A_final = zeros(n_efm, n_efm);
    A_final(final_selected_global, final_selected_global) = A_full;

    
    try
        L_final = chol(A_full + 1e-12*eye(n_sel), 'lower');
    catch
        [V, D] = eig(0.5*(A_full+A_full'));
        d = max(diag(D), 0);
        L_final = V * diag(sqrt(d)) * V';
    end

end