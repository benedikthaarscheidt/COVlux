function cluster_metrics = compute_cluster_metrics(
    clusterName, ... E_use, X_use, E_final, X_final, A_opt_QR, P_full, P_efm,
    ... reduce_colinearity, reduce_EFM_only, reduce_reactions_only,
    ... lambda_lasso, lambda_ridge, X_recon_full, L, E_red, selected_efm_idx,
    rxn_names_reduced, ... E_original_full, rxn_names_full, lambda_l21,
    max_iters, mean_influence, ... mu_final_reduced, lambda_balance)

    cluster_metrics = struct();
cluster_metrics.ClusterName = clusterName;

% -- -1. REDUCED SPACE METRICS(Data - mapped)-- -
   
        Target is E_final to count losses in the reduced space
            .[num_lost, lost_names, final_efm_idx, revived_names] =
    ... compute_unique_reaction_loss(A_opt_QR, E_final, mu_final_reduced,
                                     E_final, rxn_names_reduced);

cluster_metrics.UniqueRxnsLostCount = num_lost;

% -- -2. FULL SPACE METRICS(Projected)-- -
        % We use the SAME EFM selection(from E_final / mu_final).%
            But we target E_original_full to count losses in the full model
                .if ~isempty(E_original_full) &&
    ~isempty(rxn_names_full)

        [num_lost_full, lost_names_full, pruned_efm_idx_full,
         revived_names_full] =
    ... compute_unique_reaction_loss(A_opt_QR, E_final, mu_final_reduced,
                                     E_original_full, rxn_names_full);

% Store full space analysis data in the metrics struct instead of saving
        immediately cluster_metrics.FullSpaceAnalysis.lost_names_full =
    lost_names_full;
cluster_metrics.FullSpaceAnalysis.pruned_efm_idx_full = pruned_efm_idx_full;
cluster_metrics.FullSpaceAnalysis.revived_names_full = revived_names_full;

cluster_metrics.UniqueRxnsLostCount_Full = num_lost_full;
cluster_metrics.RevivedRxnsCount_Full = length(revived_names_full);
else cluster_metrics.UniqueRxnsLostCount_Full = NaN;
cluster_metrics.RevivedRxnsCount_Full = NaN;
end

    % 3. Compute other metrics(unchanged) cluster_metrics =
    compute_system_metrics(cluster_metrics, E_use, X_use, E_final, X_final,
                           reduce_colinearity, A_opt_QR);
cluster_metrics =
    compute_error_metrics(cluster_metrics, E_use, X_use, E_final, X_final,
                          A_opt_QR, X_recon_full, reduce_colinearity);
cluster_metrics =
    compute_fundamental_limits(cluster_metrics, E_use, X_use, E_final, X_final);
cluster_metrics = compute_solution_metrics(
    cluster_metrics, A_opt_QR, lambda_lasso, lambda_ridge, lambda_l21,
    max_iters, mean_influence, lambda_balance);
cluster_metrics = compute_projection_metrics(
    cluster_metrics, P_full, P_efm, reduce_colinearity, reduce_EFM_only,
    reduce_reactions_only);

end

    function metrics =
        compute_system_metrics(metrics, E_use, X_use, E_final, X_final,
                               reduce_colinearity, A_opt_QR)

        % SYSTEM DIMENSIONS metrics.OriginalReactions = size(E_use, 1);
metrics.OriginalEFMs = size(E_use, 2);
metrics.ReducedReactions = size(E_final, 1);
metrics.ReducedEFMs = size(E_final, 2);
metrics.ReductionRatio = size(E_final, 1) / size(E_use, 1);
metrics.EFMprunedbyL21 = sum(diag(A_opt_QR) == 0);

% ORIGINAL SYSTEM PROPERTIES metrics =
    compute_matrix_properties(metrics, X_use, 'Original');
metrics.E_full_Rank = rank(E_use);
metrics.E_full_Condition = cond(E_use);

% REDUCED SYSTEM PROPERTIES metrics =
    compute_matrix_properties(metrics, X_final, 'Reduced');
metrics.E_final_Rank = rank(E_final);
metrics.E_final_Condition = cond(E_final);

% PROBLEM STRUCTURE if reduce_colinearity metrics.Rank_Deficiency_X =
    size(X_final, 1) - rank(X_final);
metrics.Rank_Deficiency_E = size(E_final, 1) - rank(E_final);
metrics.Rank_Ratio = rank(E_final) / min(size(E_final));
metrics.DegreesOfFreedom = size(E_final, 2) ^ 2;
metrics.Equations = size(X_final, 1) ^ 2;
metrics.ProblemType = metrics.DegreesOfFreedom / metrics.Equations;
else metrics.Rank_Deficiency_X = size(X_use, 1) - rank(X_use);
metrics.Rank_Deficiency_E = size(E_use, 1) - rank(E_use);
metrics.Rank_Ratio = rank(E_use) / min(size(E_use));
metrics.DegreesOfFreedom = size(E_use, 2) ^ 2;
metrics.Equations = size(X_use, 1) ^ 2;
metrics.ProblemType = metrics.DegreesOfFreedom / metrics.Equations;
end end

    function metrics =
        compute_matrix_properties(metrics, X, prefix) %
        Compute standard matrix properties with safe error handling

            try metrics.([prefix '_X_Norm']) = norm(X, 'fro');
metrics.([prefix '_X_Trace']) = trace(X);
metrics.([prefix '_X_Det']) = det(X);

eig_vals = eig(X);
metrics.([prefix 'X_MinEig']) = min(eig_vals);
metrics.([prefix 'X_MaxEig']) = max(eig_vals);
metrics.([prefix 'X_Condition']) = cond(X);
metrics.([prefix 'X_EffectiveRank']) = sum(eig_vals > 1e-10 * max(eig_vals));

catch ME metrics.([prefix '_X_Norm']) = NaN;
metrics.([prefix '_X_Trace']) = NaN;
metrics.([prefix '_X_Det']) = NaN;
metrics.([prefix '_MinEig']) = NaN;
metrics.([prefix '_MaxEig']) = NaN;
metrics.([prefix '_Condition']) = NaN;
metrics.([prefix '_EffectiveRank']) = NaN;
    end
end


function metrics = compute_error_metrics(metrics, E_use, X_use, E_final, X_final, A_opt_QR, X_recon_full, reduce_colinearity)
    
    % --- RECONSTRUCTION SETUP ---
    % Reconstruct the covariance matrix for the reduced system
    X_recon_reduced = E_final * A_opt_QR * E_final';
    X_recon_full = E_use * A_opt_QR * E_use';
    
    % Get matrix dimensions
    m_full = size(X_use, 1);
    m_reduced = size(X_final, 1);

    % -- -1. BASELINE FROBENIUS ERRORS(Alread)-- -

        % REDUCED SYSTEM ERROR(Fit to X_final) metrics.ReducedrelativeError =
        norm(X_recon_reduced - X_final, 'fro') / norm(X_final, 'fro');
    metrics.ReducedAbsoluteError = norm(X_recon_reduced - X_final, 'fro');

    % FULL SYSTEM ERROR(Fit to X_use) metrics.FullsysRelativeError =
        norm(X_recon_full - X_use, 'fro') / norm(X_use, 'fro');
    metrics.FullsysAbsoluteError = norm(X_recon_full - X_use, 'fro');

    % Reduced System norm_1_X_final = norm(X_final, 1);
    if norm_1_X_final
      > 0 % Formula : || Residual ||
                      _F / (|| X ||
                            _1 * sqrt(m)) metrics.ReducedRelErrorMeanNorm =
          metrics.ReducedAbsoluteError / (norm_1_X_final * sqrt(m_reduced));
    else
      metrics.ReducedRelErrorMeanNorm = NaN;
    end

        % Full System norm_1_X_use = norm(X_use, 1);
    if norm_1_X_use
      > 0 metrics.FullsysRelErrorMeanNorm =
          metrics.FullsysAbsoluteError / (norm_1_X_use * sqrt(m_full));
    else
      metrics.FullsysRelErrorMeanNorm = NaN;
    end

        % -- -3. VARIANCE EXPLAINED(R ^ 2 - like Metric) 📊 -- -

       
                     % Reduced System(R ^ 2 - like)
                           norm_X_final_sq = norm(X_final, 'fro') ^ 2;
    if norm_X_final_sq
      > 1e-12 % Avoid division by zero / near - zero metrics.ReducedR2Like =
          1 - (metrics.ReducedAbsoluteError ^ 2 / norm_X_final_sq);
    else
      metrics.ReducedR2Like = 1.0;
    end

        % Full System(R ^ 2 - like) norm_X_use_sq = norm(X_use, 'fro') ^ 2;
    if norm_X_use_sq
      > 1e-12 metrics.FullsysR2Like =
          1 - (metrics.FullsysAbsoluteError ^ 2 / norm_X_use_sq);
    else
      metrics.FullsysR2Like = 1.0;
    end

        % -- -4. SPECTRAL NORM ERRORS(Worst - Case Error) 📉 -- -

        % The spectral norm is the largest singular value(norm(A, 2))

        % Reduced System norm_2_residual_reduced = norm(X_recon_reduced -
                                                            X_final,
                                                        2);
    norm_2_X_final = norm(X_final, 2);
    if norm_2_X_final
      > 0 metrics.ReducedRelErrorSpectral =
          norm_2_residual_reduced / norm_2_X_final;
    else
      metrics.ReducedRelErrorSpectral = NaN;
    end

        % Full System norm_2_residual_full = norm(X_recon_full - X_use, 2);
    norm_2_X_use = norm(X_use, 2);
    if norm_2_X_use
      > 0 metrics.FullsysRelErrorSpectral = norm_2_residual_full / norm_2_X_use;
    else
      metrics.FullsysRelErrorSpectral = NaN;
    end
    
    % --- 5. MAHALANOBIS DISTANCE RELATIVE ERROR (Advanced) 🔱 ---
    
    % Only calculate for the full system where invertibility is more likely
    % (or use the regularized/smoothed version of X_use if available)
    try
        % X must be positive definite/invertible
        X_inv_sqrt = inv(sqrtm(X_use));
    residual_weighted = X_inv_sqrt * (X_recon_full - X_use) * X_inv_sqrt;
    metrics.FullsysMahalanobisError = norm(residual_weighted, 'fro');
    catch % Failed to invert(e.g., singular or non - positive definite matrix)
                metrics.FullsysMahalanobisError = NaN;
    end

            % -- -RECONSTRUCTION ANALYSIS(Original)-- -
        if reduce_colinearity metrics.ReconstructionPenalty =
        metrics.FullsysRelativeError - metrics.ReducedrelativeError;
    metrics.ReconstructionQuality =
        metrics.ReducedrelativeError / max(metrics.FullsysRelativeError, 1e-12);
    else metrics.ReconstructionPenalty = 0;
    metrics.ReconstructionQuality = 1.0;
    end

            % -- -SCALE ANALYSIS(Original)-- -
        metrics.ReducedScaleRatio = norm(X_recon_reduced, 'fro') /
                                    norm(X_final, 'fro');
    metrics.FullScaleRatio = norm(X_recon_full, 'fro') / norm(X_use, 'fro');

    end

        function metrics =
            compute_fundamental_limits(metrics, E_use, X_use, E_final, X_final)

            % PROJECTION ERRORS metrics.ProjectionError =
                compute_projection_error(E_final, X_final);
    metrics.ProjectionError_Original = compute_projection_error(E_use, X_use);

    % SIMPLE SOLUTION ERRORS metrics.SimpleSolutionError =
        compute_simple_solution_error(E_final, X_final);
    metrics.SimpleSolutionError_Original =
        compute_simple_solution_error(E_use, X_use);

    % THEORETICAL LIMITS[metrics.TheoreticalMinError, metrics.OptimalityRatio] =
        compute_theoretical_limits(E_final, X_final,
                                   metrics.ReducedrelativeError);

    % MATRIX ALIGNMENT metrics.E_X_Alignment =
        compute_matrix_alignment(E_final, X_final);
    metrics.E_X_Alignment_Original = compute_matrix_alignment(E_use, X_use);
    end

        function projection_error =
            compute_projection_error(E, X) try[Ue, ~, ~] = svd(E, 'econ');
        X_projected = Ue * (Ue' * X * Ue) * Ue';
        projection_error = norm(X - X_projected, 'fro') / norm(X, 'fro');
    catch
        projection_error = NaN;
    end
end

function simple_error = compute_simple_solution_error(E, X)
    try
        A_simple = eye(size(E, 2));
        simple_error = norm(E * A_simple * E' - X, 'fro') / norm(X, 'fro');
    catch
        simple_error = NaN;
    end
end

function [theoretical_min, optimality_ratio] = compute_theoretical_limits(E, X, actual_error)
    try
        E_pinv = pinv(E);
        A_optimal = E_pinv * X * E_pinv';
        theoretical_min = norm(E * A_optimal * E' - X, 'fro') / norm(X, 'fro');
        optimality_ratio = actual_error / theoretical_min;
    catch
        theoretical_min = NaN;
        optimality_ratio = NaN;
    end
end

function alignment = compute_matrix_alignment(E, X)
    try
        %  compute E*E' first, then vectorize
        EET = E * E';
        vec_EET = EET(:);
        vec_X = X(:);
        alignment = (vec_EET' * vec_X) / (norm(vec_EET) * norm(vec_X));
    catch ME
        fprintf('Warning: Matrix alignment computation failed: %s\n', ME.message);
        alignment = NaN;
    end
end


function metrics = compute_solution_metrics(metrics, A_opt_QR, lambda_lasso, lambda_ridge,lambda_l21, max_iters,mean_influence,lambda_balance)

    % SOLUTION MATRIX PROPERTIES
    metrics.A_Matrix_Norm = norm(A_opt_QR, 'fro');
    metrics.A_Matrix_Trace = trace(A_opt_QR);
    
    try
        eig_A = eig(A_opt_QR);
        metrics.A_MinEig = min(eig_A);
        metrics.A_MaxEig = max(eig_A);
        metrics.A_Condition = cond(A_opt_QR);
        metrics.A_PositiveDefinite = all(eig_A > -1e-10);
    catch
        metrics.A_MinEig = NaN;
        metrics.A_MaxEig = NaN;
        metrics.A_Condition = NaN;
        metrics.A_PositiveDefinite = false;
    end
    
    
    
    % REGULARIZATION PARAMETERS
    metrics.LambdaLasso = lambda_lasso;
    metrics.LambdaRidge = lambda_ridge;
    metrics.Lambda_L21 = lambda_l21;
    metrics.Max_Iterations = max_iters;
    metrics.mean_influence = mean_influence;
    metrics.lambda_balance = lambda_balance;

end


function metrics = compute_projection_metrics(metrics, P_full, P_efm, reduce_colinearity, reduce_EFM_only, reduce_reactions_only)

    metrics.P_full_Rows = size(P_full, 1);
    metrics.P_full_Cols = size(P_full, 2);
    metrics.P_efm_Rows = size(P_efm, 1);
    metrics.P_efm_Cols = size(P_efm, 2);
    
    % PROJECTION TYPE CLASSIFICATION
    if reduce_colinearity
        if reduce_EFM_only && ~reduce_reactions_only
            metrics.ProjectionType = 'EFM_only';
        elseif reduce_reactions_only && ~reduce_EFM_only
            metrics.ProjectionType = 'Reactions_only';
        else
            metrics.ProjectionType = 'Both';
        end
    else
        metrics.ProjectionType = 'None';
    end
    
    % PROJECTION MATRIX PROPERTIES
    try
        metrics.P_full_Rank = rank(P_full);
        metrics.P_efm_Rank = rank(P_efm);
        metrics.P_full_Orthogonality = norm(P_full * P_full' - eye(size(P_full,1)), 'fro');
    catch
        metrics.P_full_Rank = NaN;
        metrics.P_efm_Rank = NaN;
        metrics.P_full_Orthogonality = NaN;
    end
end
