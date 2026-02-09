function [A_psd, L_final] = solve_psd_final_qr_psd(E_sub, X,verbose, saveDir, clusterName)
    % Robust PSD solution with guaranteed Cholesky factorization
    % Inputs:
    %   verbose:     (Optional) boolean to toggle debug prints. Default is false.
    %   saveDir:     (Optional) Directory to save the histogram.
    %   clusterName: (Optional) Name prefix for the plot file.
    
    if nargin < 3, verbose = false; end
    if nargin < 4, saveDir = ''; end
    if nargin < 5, clusterName = 'UnknownCluster'; end

    % Get optimal solution
    A_opt = solve_QR_regularized_new(E_sub, X, 1e-3);
    
    % Try Cholesky first
    [L, status] = chol(A_opt, 'lower');
    
    if status == 0
        % Matrix is already PSD
        A_psd = A_opt;
        if verbose
            fprintf('PSD QR: Matrix already PSD, Cholesky_status=%d\n', status);
        end
    else
        if verbose
            fprintf('WARNING: Cholesky failed initially (status=%d)\n', status);
        end
        
        % Use robust PSD projection
        % Note: Passing 1e-12 as default ratio to match signature if needed
        A_psd = robust_psd_projection(A_opt, 1e-12, verbose);
        
        % Verify Cholesky works on the projected matrix
        [L, final_status] = chol(A_psd, 'lower');
        
        if final_status > 0
            if verbose
                fprintf('EMERGENCY: Cholesky still failing, using regularized identity\n');
            end
            % Last resort - use regularized identity
            A_psd = eye(size(A_opt)) * norm(X, 'fro') / norm(E_sub * E_sub', 'fro');
            [L, final_status] = chol(A_psd, 'lower');
        end
        
        if verbose
            fprintf('PSD QR: After projection, Cholesky_status=%d\n', final_status);
        end
    end
    
    % Final validation
    [L_final, status_final] = chol(A_psd, 'lower');
    if status_final > 0
        error('FATAL: Matrix still not PSD after all corrections');
    end
    
    if verbose
        error_val = norm(E_sub * A_psd * E_sub' - X, 'fro') / norm(X, 'fro');
        fprintf('After PSD mapping QR: Final error=%.6f, Cholesky_status=%d\n', error_val, status_final);
    end

    % --- NEW: HISTOGRAM PLOTTING ---
    if ~isempty(saveDir)
        try
            fig = figure('Name', 'L_Matrix_Distribution', 'Color', 'w', 'Visible', 'off');
            
            % Extract non-zero values (L is lower triangular)
            l_vals = L_final(L_final ~= 0);
            
            histogram(l_vals, 50, 'FaceColor', [0.2 0.6 0.5], 'EdgeColor', 'none');
            grid on;
            title(['Distribution of L Entries: ' clusterName], 'Interpreter', 'none');
            xlabel('Value');
            ylabel('Frequency');
            subtitle(sprintf('Min: %.2e | Max: %.2e', min(l_vals), max(l_vals)));
            
            % Generate filename
            safeName = regexprep(clusterName, '[^a-zA-Z0-9]', '_');
            fName = fullfile(saveDir, sprintf('%s_L_distribution.png', safeName));
            
            saveas(fig, fName);
            if verbose, fprintf('  Saved L-matrix histogram to: %s\n', fName); end
            close(fig);
        catch ME
            fprintf('  Warning: Could not save L histogram: %s\n', ME.message);
            if exist('fig', 'var'), close(fig); end
        end
    end
end