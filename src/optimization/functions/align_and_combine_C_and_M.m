function [X_combined, rxns_common, stats, lambda_balance] = align_and_combine_C_and_M(T_cov, T_outer, saveDir, clusterName, verbose)
% ALIGN_AND_COMBINE_C_AND_M
% Inputs:
%   T_cov:       Table containing the Covariance Matrix (C)
%   T_outer:     Table containing the Mean Outer Product (M = mu*mu'), or []
%   saveDir:     (Optional) Directory to save the plot
%   clusterName: (Optional) Name of the cluster for file naming
%   verbose:     (Optional) Boolean to toggle debug prints. Default is false.
%
% Outputs:
%   X_combined:  Matrix X = C + lambda * M
%   rxns_common: String array of reaction IDs present in X_combined
%   stats:       Struct containing lambda and energies for logging

    if nargin < 5, verbose = false; end
    if nargin < 3, saveDir = ''; end
    if nargin < 4, clusterName = 'UnknownCluster'; end

    % 1. Extract Covariance Data
    rxn_C = string(T_cov.Properties.RowNames);
    C_matrix = table2array(T_cov);
    
    % 2. Check if Outer Product exists
    if isempty(T_outer)
        % Fallback: Covariance Only
        X_combined = C_matrix;
        rxns_common = rxn_C;
        
        lambda_balance = 0; 
        
        stats.lambda = 0;
        stats.energy_C = norm(C_matrix, 'fro');
        stats.energy_M = 0;
        mode_str = 'Covariance Only';
    else
        % 3. Extract Outer Product Data
        rxn_M = string(T_outer.Properties.RowNames);
        M_matrix = table2array(T_outer);
        
        % 4. ALIGNMENT (Intersection)
        [rxns_common, idx_C, idx_M] = intersect(rxn_C, rxn_M, 'stable');
        
        if verbose && (length(rxns_common) ~= length(rxn_C))
            fprintf('  [Align] Covariance Rxns: %d | Mean Rxns: %d -> Common: %d\n', ...
                length(rxn_C), length(rxn_M), length(rxns_common));
        end
        
        C_aligned = C_matrix(idx_C, idx_C);
        M_aligned = M_matrix(idx_M, idx_M);
        
        % 5. THE EQUALIZER (Calculate Lambda)
        energy_C = norm(C_aligned, 'fro');
        energy_M = norm(M_aligned, 'fro');
        
        % (Your logic override was technically here, keeping it as written)
        if energy_M > 1e-12
            lambda_balance= energy_C / energy_M;
        else
            lambda_balance = 0; 
        end
        
        
        %lambda_balance = 1; 
        
        if verbose
            fprintf('  [Balancing] Energy C: %.2e | Energy M: %.2e\n', energy_C, energy_M);
            fprintf('  [Balancing] Applied Weight (lambda): %.2e\n', lambda_balance);
        end
        
        % 6. COMBINE (Second Moment)
        X_combined = C_aligned + (lambda_balance * M_aligned);
        

        fprintf('Check for Second moment Negativity: Smallest Value = %.4e | Count < 0 = %d\n', min(X_combined(:)), sum(X_combined(:) < 0));
        %X_combined(X_combined<0)=0;
        stats.lambda = lambda_balance;
        stats.energy_C = energy_C;
        stats.energy_M = energy_M;
        mode_str = sprintf('Second Moment (Lambda=%.2e)', lambda_balance);
    end
    
    % --- 7. PLOTTING & SAVING ---
    if verbose && ~isempty(saveDir)
        try
            fig = figure('Name', 'Target Matrix Distribution', 'Color', 'w', 'Visible', 'off');
            
            % Plot histogram (using log scale for Y to see the tail)
            histogram(X_combined(:), 100, 'FaceColor', [0.6 0.2 0.6], 'EdgeColor', 'none');
            set(gca, 'YScale', 'log');
            
            title({'Distribution of Second Moment Matrix Entries', '(Balanced Target X)'});
            xlabel('Value');
            ylabel('Count (Log Scale)');
            grid on;
            
            subtitle(sprintf('Mode: %s | Max Val: %.2e', mode_str, max(abs(X_combined(:)))));
            
            filename = fullfile(saveDir, sprintf('%s_second_moment_distribution.png', clusterName));
            
            saveas(fig, filename);
            fprintf('  Saved second moment distribution plot to: %s\n', filename);
            close(fig);
        catch ME
            fprintf('  Warning: Could not save second moment distribution: %s\n', ME.message);
            if exist('fig', 'var'), close(fig); end
        end
    end
end