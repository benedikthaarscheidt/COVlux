function [A_opt, E, L, metrics] = closed_form_A(E, X, rxnNames, verbose, plotDir, clusterName, lambda_l21, div_by_reactions, threshold_method)
    % Default to diagonal if threshold_method is not provided in the function call
    if nargin < 9
        threshold_method = 'diagonal'; 
    end

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

    sym_discrepancy = norm(A - A', 'fro');


    eig_vals = eig(A);
    min_eig = min(eig_vals);
    fprintf('Minimum Eigenvalue:   %e\n', min_eig);
    if sym_discrepancy < 1e-12
        fprintf('Result: A is numerically SYMMETRIC. ✓\n');
    else
        fprintf('Result: A is ASYMMETRIC (Error: %e). ✗\n', sym_discrepancy);
        A = (A + A') * 0.5;
    end
    
    
    f3 = figure('Visible', 'off');
    histogram(diag(A), 50);
    title(sprintf('Matrix A Values (%s)', clusterName), 'Interpreter', 'none');
    saveas(f3, fullfile(plotDir, sprintf('%s_Matrix_A.png', clusterName)));
    close(f3);
    
    % --- NEW: Toggleable Metric Calculation ---
    if strcmpi(threshold_method, 'rowsum')
        fprintf("Metric: Absolute Off-Diagonal Row Sums \n");
        % Subtract the diagonal from A, then take the absolute sum of each row
        A_no_diag = A - diag(diag(A));
        metric_values = sum(abs(A_no_diag), 2);
        metric_name = 'RowSums';
    else
        fprintf("Metric: Matrix Diagonal (Variances) \n");
        % Use the standard intrinsic variance
        metric_values = diag(A);
        metric_name = 'Diagonal';
    end
    
    % --- Apply Dynamic Thresholding based on chosen metric ---
    max_val = max(metric_values);
    dynamic_threshold = max_val * lambda_l21; 
    
    fprintf("Max %s is %f. Dynamic threshold set to: %f \n", metric_name, max_val, dynamic_threshold);
    
    valid_mask = metric_values > 0;
    surviving_values = metric_values(valid_mask);
    
    f4 = figure('Visible', 'off');
    if ~isempty(surviving_values)
        % Transform the values to log10 scale
        log_surviving_values = log10(surviving_values);
        
        % Generate the density plot (Probability Density Function)
        histogram(log_surviving_values, 50, 'Normalization', 'pdf');
        hold on; % Keep the histogram while we draw the line
        
        % Calculate the log of the threshold
        log_thresh = log10(dynamic_threshold);
        
        % Add the vertical line at the threshold
        % 'LineWidth', 2 makes it stand out against the bars
        xl = xline(log_thresh, '--r', 'Threshold', 'LineWidth', 2, 'LabelOrientation', 'aligned');
        xl.FontSize = 10;
        
        % Update title and labels
        title_str = sprintf('Log10 Density Dist. (%s): Threshold at %.2e', clusterName, dynamic_threshold);
        title(title_str, 'Interpreter', 'none');
        xlabel(['log10(', metric_name, ')']);
        ylabel('Probability Density');
        
        hold off;
    else
        title(sprintf('No %s > 0 found for (%s)', metric_name, clusterName), 'Interpreter', 'none');
        text(0.5, 0.5, 'All EFMs are zero', 'HorizontalAlignment', 'center');
    end
    
    saveas(f4, fullfile(plotDir, sprintf('%s_Matrix_A_%s_log_density_with_threshold.png', clusterName, metric_name)));
    close(f4);
    
    selected_efm_idx = find(metric_values > dynamic_threshold);
    fprintf("Number of selected EFMs: %d \n", length(selected_efm_idx));
    
    A_pruned = A(selected_efm_idx, selected_efm_idx);
    
    A_opt = zeros(size(A)); 
    A_opt(selected_efm_idx, selected_efm_idx) = A_pruned;
    
    fprintf("Smallest value in diag(A): %f \n", min(diag(A)));
    fprintf("Largest value in diag(A): %f \n", max(diag(A)));
    fprintf("Smallest value in A: %f \n", min(A(:)));
    fprintf("Largest value in A: %f \n", max(A(:)));
    
    fprintf("########### after thresholding A ########### \n");
    if ~isempty(A_pruned)
        fprintf("Smallest value in diag(A_pruned): %f \n", min(diag(A_pruned)));
        fprintf("Largest value in diag(A_pruned): %f \n", max(diag(A_pruned)));
        fprintf("Smallest value in A_pruned: %f \n", min(A_pruned(:)));
        fprintf("Largest value in A_pruned: %f \n", max(A_pruned(:)));
    else
        fprintf("A_pruned is empty (all EFMs removed).\n");
    end
    
    metrics = struct();
    
    A_ridge = A_opt + eye(size(A_opt)) * 1e-8;
    L = chol(A_ridge, 'lower');
    
    f5 = figure('Visible', 'off');
    histogram(L(:), 50);
    title(sprintf('Matrix L Values (%s)', clusterName), 'Interpreter', 'none');
    saveas(f5, fullfile(plotDir, sprintf('%s_Matrix_L.png', clusterName)));
    close(f5);
end