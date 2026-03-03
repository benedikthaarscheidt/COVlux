function [A_opt, E, L, metrics] = closed_form_A_cov(E, X, rxnNames, verbose, plotDir, clusterName, BETA, SPARSITY_WEIGHT, div_by_reactions)
    fprintf('\n=== Calculate A via closed form ===\n');
    if ~exist(plotDir, 'dir')
        mkdir(plotDir);
    end

    EtE = E' * E;
    [r, c] = size(EtE);
    fprintf("Dim of EtE: %d x %d \n", r, c);
    [U, S, V] = svd(EtE, 'econ');
    s_vals = diag(S);
    f1 = figure('Visible', 'off');
    histogram(s_vals, 50);         
    title(sprintf('Singular Values (%s)', clusterName), 'Interpreter', 'none');
    saveas(f1, fullfile(plotDir, sprintf('%s_singular_values.png', clusterName)));
    close(f1); 
    max_s = max(s_vals);
    fprintf("Max singular value: %f \n", max_s);
    fprintf("Min singular value: %f \n", min(s_vals));
    tolerance = 1e-6 * max_s; 
    keep_idx = s_vals > tolerance;
    fprintf("Number of selected singular values: %d \n", sum(keep_idx));
    s_inv = zeros(size(s_vals));
    s_inv(keep_idx) = 1 ./ s_vals(keep_idx);
    EtE_pseudoinv = V * diag(s_inv) * U';
    f2 = figure('Visible', 'off');
    histogram(EtE_pseudoinv(:), 50);
    title(sprintf('EtE Pseudoinverse Values (%s)', clusterName), 'Interpreter', 'none');
    saveas(f2, fullfile(plotDir, sprintf('%s_EtE_pseudoinv.png', clusterName)));
    close(f2);
    A = EtE_pseudoinv * E' * X * E * EtE_pseudoinv;
    A = (A + A') * 0.5;
    f3 = figure('Visible', 'off');
    histogram(diag(A), 50);
    title(sprintf('Matrix A Values (%s)', clusterName), 'Interpreter', 'none');
    saveas(f3, fullfile(plotDir, sprintf('%s_Matrix_A.png', clusterName)));
    close(f3);


    variances = diag(A);
    max_var = max(variances);
    dynamic_threshold = max_var * SPARSITY_WEIGHT; 
    
    fprintf("Max variance is %f. Dynamic threshold set to: %f \n", max_var, dynamic_threshold);
    
    surviving_variances = variances(variances > dynamic_threshold);
    
    f4 = figure('Visible', 'off');
    if ~isempty(surviving_variances)
        histogram(surviving_variances, 50);
        title(sprintf('Matrix A Diagonal > %.4f (%s)', dynamic_threshold, clusterName), 'Interpreter', 'none');
    else
        title(sprintf('No variances > %.4f (%s)', dynamic_threshold, clusterName), 'Interpreter', 'none');
        text(0.5, 0.5, 'All EFMs pruned by dynamic cut', 'HorizontalAlignment', 'center');
    end
    saveas(f4, fullfile(plotDir, sprintf('%s_Matrix_A_diag_cut.png', clusterName)));
    close(f4);
    
    % Apply the exact same dynamic threshold to isolate the EFMs
    selected_efm_idx = find(variances > dynamic_threshold);
    fprintf("Number of selected EFMs: %d \n", length(selected_efm_idx));
    
    A_pruned = A(selected_efm_idx, selected_efm_idx);
    
    A_opt = zeros(size(A)); % zeros(size(A)) is safer than zeros(size(A,2))
    A_opt(selected_efm_idx, selected_efm_idx) = A_pruned;
    
    fprintf("Smallest value in diag(A): %f \n", min(diag(A)));
    fprintf("Largest value in diag(A): %f \n", max(diag(A)));
    fprintf("Smallest value in A: %f \n", min(A(:)));
    fprintf("Largest value in A: %f \n", max(A(:)));
    
    fprintf("########### after thresholding A ########### \n");
    fprintf("Smallest value in diag(A_pruned): %f \n", min(diag(A_pruned)));
    fprintf("Largest value in diag(A_pruned): %f \n", max(diag(A_pruned)));
    fprintf("Smallest value in A_pruned: %f \n", min(A_pruned(:)));
    fprintf("Largest value in A_pruned: %f \n", max(A_pruned(:)));
    
    metrics = struct();
    
    % --- CRITICAL FIX: Changed A to A_opt for the Cholesky ridge ---
    A_ridge = A_opt + eye(size(A_opt)) * 1e-8;
    L = chol(A_ridge, 'lower');
    
    f5 = figure('Visible', 'off');
    histogram(L(:), 50);
    title(sprintf('Matrix L Values (%s)', clusterName), 'Interpreter', 'none');
    saveas(f5, fullfile(plotDir, sprintf('%s_Matrix_L.png', clusterName)));
    close(f5);
end