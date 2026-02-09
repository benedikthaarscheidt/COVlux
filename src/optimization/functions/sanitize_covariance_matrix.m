function [X_clean, stats] = sanitize_covariance_matrix(X_raw, threshold, verbose, saveDir, clusterName)
% SANITIZE_COVARIANCE_MATRIX with saving capability
    
    if nargin < 3, verbose = true; end
    if nargin < 2, threshold = []; end 
    
    X_clean = X_raw;
    n = size(X_clean, 1);
    
    % --- Safety Check for NaNs/Infs ---
    if any(isnan(X_clean(:))) || any(isinf(X_clean(:)))
        if verbose
            fprintf('  ⚠️ [Safety] Input matrix contains NaNs or Infs. Replacing with 0.\n');
        end
        X_clean(isnan(X_clean)) = 0;
        X_clean(isinf(X_clean)) = 0;
    end
    
    % --- 1. Identify Off-Diagonal Elements ---
    off_diag_mask = ~eye(n);
    off_diag_vals = abs(X_clean(off_diag_mask));
    
    % --- 2. Determine Threshold ---
    if isempty(threshold)
        % Auto Mode: 25th percentile of NON-ZERO values
        non_zero_vals = off_diag_vals(off_diag_vals > 0);
        if isempty(non_zero_vals)
            actual_cutoff = 0;
        else
            actual_cutoff = prctile(non_zero_vals, 25);
        end
        mode_str = 'Auto (25% quantile)';
    else
        % Manual Mode
        actual_cutoff = threshold;
        mode_str = 'Manual';
    end
    
    % --- 3. Apply Cleanup (Protecting Diagonal) ---
    mask_noise = abs(X_clean) < actual_cutoff;
    
    % Protect Diagonal (Variance)
    diag_indices = 1:(n+1):numel(X_clean);
    mask_noise(diag_indices) = false;
    
    X_clean(mask_noise) = 0;
    
    n_removed = sum(mask_noise(:));
    
    % --- 4. Stability (Ridge) ---
    X_clean = X_clean + 1e-9 * eye(n);
    X_clean = project_to_psd(X_clean, 1e-8);
    
    % --- 5. Logging & Visualization ---
    if verbose
        fprintf('\n[Target Covariance X Sanitization]\n');
        fprintf('  Mode            : %s\n', mode_str);
        fprintf('  Cutoff Value    : %.2e\n', actual_cutoff);
        fprintf('  Noise Removed   : %d elements (%.1f%% of matrix)\n', ...
            n_removed, (n_removed/numel(X_clean))*100);
        fprintf('  ------------------------------------------\n');
        
        % --- VISUALIZATION BLOCK ---
        fig = figure('Name', 'Covariance Distribution Analysis', 'Color', 'w', 'Visible', 'off');
        
        % Filter out perfect zeros for log-plot clarity
        plot_vals = off_diag_vals(off_diag_vals > 1e-20);
        
        if isempty(plot_vals)
            title('Matrix is empty or all zeros');
        else
            % Plot Histogram of Log10 magnitudes
            histogram(log10(plot_vals), 100, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
            hold on;
            
            % Draw the Threshold Line
            xline(log10(actual_cutoff), '--r', 'LineWidth', 2, ...
                'Label', sprintf('Cutoff: 10^{%.1f}', log10(actual_cutoff)));
            
            % Formatting
            grid on;
            xlabel('Log_{10} |Covariance|');
            ylabel('Count (Frequency)');
            title({'Distribution of Off-Diagonal Covariances', ...
                   ['(Values < 10^{', sprintf('%.1f', log10(actual_cutoff)), '} removed)']});
            
            subtitle(sprintf('Mode: %s | Removed: %.1f%%', mode_str, (n_removed/numel(X_clean))*100));
            legend('Data Distribution', 'Sanitization Threshold');
            hold off;
        end
        
        % Save figure if directories are provided
        if nargin >= 4 && ~isempty(saveDir)
            if nargin >= 5
                filename = fullfile(saveDir, sprintf('%s_covariance_distribution.png', clusterName));
            else
                filename = fullfile(saveDir, 'covariance_distribution.png');
            end
            
            % Save in multiple formats
            saveas(fig, filename);
            saveas(fig, strrep(filename, '.png', '.pdf'));
            saveas(fig, strrep(filename, '.png', '.fig'));
            
            fprintf('Saved covariance distribution to: %s\n', filename);
        end
        
        close(fig);
    end
    
    stats.cutoff_used = actual_cutoff;
    stats.removed = n_removed;
end