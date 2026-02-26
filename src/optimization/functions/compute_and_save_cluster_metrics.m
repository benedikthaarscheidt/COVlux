function cluster_metrics = compute_and_save_cluster_metrics(resultsDir, clusterName, ...
    E_use, X_use, E_final, X_final, A_opt_QR, P_full, P_efm, ...
    reduce_colinearity, reduce_EFM_only, reduce_reactions_only, ...
    lambda_lasso, lambda_ridge, X_recon_full, L, E_red, selected_efm_idx, rxn_names_reduced, ...
    E_original_full, rxn_names_full, lambda_l21, max_iters, mean_influence, ...
    mu_final_reduced,lambda_balance)  
    
    cluster_metrics = struct();
    cluster_metrics.ClusterName = clusterName;
    
    
    [num_lost, lost_names, final_efm_idx, revived_names] = ...
        compute_unique_reaction_loss(A_opt_QR, E_final, mu_final_reduced, E_final, rxn_names_reduced,false);
    
    cluster_metrics.UniqueRxnsLostCount = num_lost;
    
    
    if ~isempty(E_original_full) && ~isempty(rxn_names_full)
        saveDir=fullfile(resultsDir,"plots");
        [num_lost_full, lost_names_full, pruned_efm_idx_full, revived_names_full] = ...
            compute_unique_reaction_loss(A_opt_QR, E_final, mu_final_reduced, E_original_full, rxn_names_full,true);

        
        
        % Save full space files
        save_full_space_analysis(resultsDir, clusterName, ...
            lost_names_full, pruned_efm_idx_full, selected_efm_idx, revived_names_full);
        
        cluster_metrics.UniqueRxnsLostCount_Full = num_lost_full;
        cluster_metrics.RevivedRxnsCount_Full = length(revived_names_full);
    else
        cluster_metrics.UniqueRxnsLostCount_Full = NaN;
        cluster_metrics.RevivedRxnsCount_Full = NaN;
    end
    
    % 3. Compute other metrics (unchanged)
    cluster_metrics = compute_system_metrics(cluster_metrics, E_use, X_use, E_final, X_final, reduce_colinearity, A_opt_QR);
    cluster_metrics = compute_error_metrics(cluster_metrics, E_use, X_use, E_final, X_final, A_opt_QR, X_recon_full, reduce_colinearity);
    cluster_metrics = compute_fundamental_limits(cluster_metrics, E_use, X_use, E_final, X_final);
    cluster_metrics = compute_solution_metrics(cluster_metrics, A_opt_QR, lambda_lasso, lambda_ridge, lambda_l21, max_iters, mean_influence,lambda_balance);
    cluster_metrics = compute_projection_metrics(cluster_metrics, P_full, P_efm, reduce_colinearity, reduce_EFM_only, reduce_reactions_only);
    
    % 4. Save results
    save_cluster_results(resultsDir, clusterName, cluster_metrics, A_opt_QR, P_full, P_efm, ...
        E_use, X_use, E_final, X_final, L, E_red, selected_efm_idx, lost_names_full);
    
    
end

function metrics = compute_system_metrics(metrics, E_use, X_use, E_final, X_final, reduce_colinearity,A_opt_QR)

    % SYSTEM DIMENSIONS
    metrics.OriginalReactions = size(E_use, 1);
    metrics.OriginalEFMs = size(E_use, 2);
    metrics.ReducedReactions = size(E_final, 1);
    metrics.ReducedEFMs = size(E_final, 2);
    metrics.ReductionRatio = size(E_final, 1) / size(E_use, 1);
    metrics.EFMprunedbyL21 = sum(diag(A_opt_QR) ==0);
    
     
    
    % ORIGINAL SYSTEM PROPERTIES
    metrics = compute_matrix_properties(metrics, X_use, 'Original');
    metrics.E_full_Rank = rank(E_use);
    metrics.E_full_Condition = cond(E_use);
    
    % REDUCED SYSTEM PROPERTIES  
    metrics = compute_matrix_properties(metrics, X_final, 'Reduced');
    metrics.E_final_Rank = rank(E_final);
    metrics.E_final_Condition = cond(E_final);
    
    % PROBLEM STRUCTURE
    if reduce_colinearity
        metrics.Rank_Deficiency_X = size(X_final, 1) - rank(X_final);
        metrics.Rank_Deficiency_E = size(E_final, 1) - rank(E_final);
        metrics.Rank_Ratio = rank(E_final) / min(size(E_final));
        metrics.DegreesOfFreedom = size(E_final, 2)^2;
        metrics.Equations = size(X_final, 1)^2;
        metrics.ProblemType = metrics.DegreesOfFreedom / metrics.Equations;
    else
        metrics.Rank_Deficiency_X = size(X_use, 1) - rank(X_use);
        metrics.Rank_Deficiency_E = size(E_use, 1) - rank(E_use);
        metrics.Rank_Ratio = rank(E_use) / min(size(E_use));
        metrics.DegreesOfFreedom = size(E_use, 2)^2;
        metrics.Equations = size(X_use, 1)^2;
        metrics.ProblemType = metrics.DegreesOfFreedom / metrics.Equations;
    end
end


function metrics = compute_matrix_properties(metrics, X, prefix)
% Compute standard matrix properties with safe error handling

    try
        metrics.([prefix '_X_Norm']) = norm(X, 'fro');
        metrics.([prefix '_X_Trace']) = trace(X);
        metrics.([prefix '_X_Det']) = det(X);
        
        eig_vals = eig(X);
        metrics.([prefix 'X_MinEig']) = min(eig_vals);
        metrics.([prefix 'X_MaxEig']) = max(eig_vals);
        metrics.([prefix 'X_Condition']) = cond(X);
        metrics.([prefix 'X_EffectiveRank']) = sum(eig_vals > 1e-10 * max(eig_vals));
        
    catch ME
        metrics.([prefix '_X_Norm']) = NaN;
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
    
    % --- 1. BASELINE FROBENIUS ERRORS ---
    
    % REDUCED SYSTEM ERROR (Fit to X_final)
    metrics.ReducedrelativeError = norm(X_recon_reduced - X_final, 'fro') / norm(X_final, 'fro');
    metrics.ReducedAbsoluteError = norm(X_recon_reduced - X_final, 'fro');
    
    % FULL SYSTEM ERROR (Fit to X_use)
    metrics.FullsysRelativeError = norm(X_recon_full - X_use, 'fro') / norm(X_use, 'fro');
    metrics.FullsysAbsoluteError = norm(X_recon_full - X_use, 'fro');
    
   
    % Reduced System
    norm_1_X_final = norm(X_final, 1);
    if norm_1_X_final > 0
        % Formula: ||Residual||_F / (||X||_1 * sqrt(m))
        metrics.ReducedRelErrorMeanNorm = metrics.ReducedAbsoluteError / (norm_1_X_final * sqrt(m_reduced));
    else
        metrics.ReducedRelErrorMeanNorm = NaN;
    end
    
    % Full System
    norm_1_X_use = norm(X_use, 1);
    if norm_1_X_use > 0
        metrics.FullsysRelErrorMeanNorm = metrics.FullsysAbsoluteError / (norm_1_X_use * sqrt(m_full));
    else
        metrics.FullsysRelErrorMeanNorm = NaN;
    end
  
    
    % Reduced System (R^2-like)
    norm_X_final_sq = norm(X_final, 'fro')^2;
    if norm_X_final_sq > 1e-12 % Avoid division by zero/near-zero
        metrics.ReducedR2Like = 1 - (metrics.ReducedAbsoluteError^2 / norm_X_final_sq);
    else
        metrics.ReducedR2Like = 1.0;
    end
    
    % Full System (R^2-like)
    norm_X_use_sq = norm(X_use, 'fro')^2;
    if norm_X_use_sq > 1e-12
        metrics.FullsysR2Like = 1 - (metrics.FullsysAbsoluteError^2 / norm_X_use_sq);
    else
        metrics.FullsysR2Like = 1.0;
    end
    
    % --- 4. SPECTRAL NORM ERRORS (Worst-Case Error) 📉 ---
    
    % The spectral norm is the largest singular value (norm(A, 2))
    
    % Reduced System
    norm_2_residual_reduced = norm(X_recon_reduced - X_final, 2);
    norm_2_X_final = norm(X_final, 2);
    if norm_2_X_final > 0
        metrics.ReducedRelErrorSpectral = norm_2_residual_reduced / norm_2_X_final;
    else
        metrics.ReducedRelErrorSpectral = NaN;
    end
    
    % Full System
    norm_2_residual_full = norm(X_recon_full - X_use, 2);
    norm_2_X_use = norm(X_use, 2);
    if norm_2_X_use > 0
        metrics.FullsysRelErrorSpectral = norm_2_residual_full / norm_2_X_use;
    else
        metrics.FullsysRelErrorSpectral = NaN;
    end
    
    
    try
       
        X_inv_sqrt = inv(sqrtm(X_use)); 
        residual_weighted = X_inv_sqrt * (X_recon_full - X_use) * X_inv_sqrt;
        metrics.FullsysMahalanobisError = norm(residual_weighted, 'fro');
    catch
        % Failed to invert (e.g., singular or non-positive definite matrix)
        metrics.FullsysMahalanobisError = NaN;
    end
    
    % --- RECONSTRUCTION ANALYSIS ---
    if reduce_colinearity
        metrics.ReconstructionPenalty = metrics.FullsysRelativeError - metrics.ReducedrelativeError;
        metrics.ReconstructionQuality = metrics.ReducedrelativeError / max(metrics.FullsysRelativeError, 1e-12);
    else
        metrics.ReconstructionPenalty = 0;
        metrics.ReconstructionQuality = 1.0;
    end
    
    % --- SCALE ANALYSIS ---
    metrics.ReducedScaleRatio = norm(X_recon_reduced, 'fro') / norm(X_final, 'fro');
    metrics.FullScaleRatio = norm(X_recon_full, 'fro') / norm(X_use, 'fro');
    
end



function metrics = compute_fundamental_limits(metrics, E_use, X_use, E_final, X_final)

    % PROJECTION ERRORS
    metrics.ProjectionError = compute_projection_error(E_final, X_final);
    metrics.ProjectionError_Original = compute_projection_error(E_use, X_use);
    
    % SIMPLE SOLUTION ERRORS
    metrics.SimpleSolutionError = compute_simple_solution_error(E_final, X_final);
    metrics.SimpleSolutionError_Original = compute_simple_solution_error(E_use, X_use);
    
    % THEORETICAL LIMITS
    [metrics.TheoreticalMinError, metrics.OptimalityRatio] = compute_theoretical_limits(E_final, X_final, metrics.ReducedrelativeError);
    
    % MATRIX ALIGNMENT
    metrics.E_X_Alignment = compute_matrix_alignment(E_final, X_final);
    metrics.E_X_Alignment_Original = compute_matrix_alignment(E_use, X_use);
end

function projection_error = compute_projection_error(E, X)
    try
        [Ue, ~, ~] = svd(E, 'econ');
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
        % Correct syntax - compute E*E' first, then vectorize
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

function save_cluster_results(resultsDir, clusterName, metrics, A_opt_QR, P_full, P_efm, ...
    E_use, X_use, E_final, X_final, L, E_red, selected_efm_idx, lost_names)
% SAVE_CLUSTER_RESULTS - Saves all main matrices and metrics for a cluster run.
% Added saving of unique lost reaction names.
    try
        % Ensure directory exists
        if ~exist(resultsDir, 'dir')
            mkdir(resultsDir);
        end
        
        % 1. Save metrics as CSV
        metrics_table = struct2table(metrics);
        metrics_file = fullfile(resultsDir, sprintf('%s_metrics.csv', clusterName));
        writetable(metrics_table, metrics_file);
        
        % 2. Save main matrices
        outFileA = fullfile(resultsDir, sprintf('%s_A.csv', clusterName));
        outFileP = fullfile(resultsDir, sprintf('%s_P_rxn.csv', clusterName));
        outFileL = fullfile(resultsDir, sprintf('%s_L.csv', clusterName));
        outfileEfinal = fullfile(resultsDir, sprintf('%s_E_final.csv', clusterName));
        outfileEred = fullfile(resultsDir, sprintf('%s_E_red.csv', clusterName));
        outfileSelectedEFM = fullfile(resultsDir, sprintf('%s_selected_efm_idx.csv', clusterName));
        outFileLostRxns = fullfile(resultsDir, sprintf('%s_unique_lost_rxns.csv', clusterName)); % NEW
        
        if exist('writematrix', 'file')
            writematrix(A_opt_QR, outFileA);
            writematrix(P_full, outFileP);
            writematrix(L, outFileL);
            writematrix(E_final, outfileEfinal);
            writematrix(E_red, outfileEred);
            writematrix(selected_efm_idx, outfileSelectedEFM);
            writetable(table(lost_names), outFileLostRxns); % Save lost names as a table/CSV
            
            if ~isempty(P_efm) && ~isequal(P_efm, eye(size(P_efm)))
                outFileP_efm = fullfile(resultsDir, sprintf('%s_P_efm.csv', clusterName));
                writematrix(P_efm, outFileP_efm);
            end
        else
            csvwrite(outFileA, A_opt_QR);
            csvwrite(outFileP, P_full);
            csvwrite(outFileL, L);
            csvwrite(outfileEfinal, E_final);
            csvwrite(outfileEred, E_red);
            csvwrite(outfileSelectedEFM, selected_efm_idx);
            csvwrite(outFileLostRxns, lost_names); % Fallback for older MATLAB versions
            
            if ~isempty(P_efm) && ~isequal(P_efm, eye(size(P_efm)))
                outFileP_efm = fullfile(resultsDir, sprintf('%s_P_efm.csv', clusterName));
                csvwrite(outFileP_efm, P_efm);
            end
        end
        
        % 3. Save debug matrices
        debug_file = fullfile(resultsDir, sprintf('%s_debug.mat', clusterName));
        save(debug_file, 'E_use', 'X_use', 'E_final', 'X_final', 'A_opt_QR', 'L', ...
            'P_full', 'P_efm', 'E_red', 'selected_efm_idx', 'metrics', 'lost_names');
        
        fprintf('✓ Saved comprehensive results for %s\n', clusterName);
        
    catch ME
        fprintf('✗ ERROR saving results for %s: %s\n', clusterName, ME.message);
    end
end





function save_full_space_analysis(resultsDir, clusterName, ...
    lost_names_full, pruned_efm_idx_full, selected_efm_idx_full, revived_names_full)
    
    try
        % 1. Save full space lost reactions
        outFileLostRxnsFull = fullfile(resultsDir, sprintf('%s_unique_lost_rxns_full.csv', clusterName));
        writetable(table(lost_names_full), outFileLostRxnsFull);
        
        % 2. Save full space EFM indices
        outFilePrunedEFMFull = fullfile(resultsDir, sprintf('%s_pruned_efm_idx_full.csv', clusterName));
        outFileSelectedEFMFull = fullfile(resultsDir, sprintf('%s_selected_efm_idx_full.csv', clusterName));
        
        writematrix(pruned_efm_idx_full, outFilePrunedEFMFull);
        writematrix(selected_efm_idx_full, outFileSelectedEFMFull);
        
        % 3. Save Revived Reactions (NEW)
        if nargin >= 6 && ~isempty(revived_names_full)
            outFileRevivedFull = fullfile(resultsDir, sprintf('%s_revived_rxns_full.csv', clusterName));
            writetable(table(revived_names_full), outFileRevivedFull);
            fprintf('  Saved revived reactions list for %s\n', clusterName);
        end
        
        fprintf('  Saved full space analysis files for %s\n', clusterName);
    catch ME
        fprintf('  Warning: Could not save full space analysis: %s\n', ME.message);
    end
end