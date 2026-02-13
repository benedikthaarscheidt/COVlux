function [X_clean, scale_factor] = prepare_covariance_matrix(X, target_max, verbose)
    % Prepare X for optimization: Clip -> Scale -> Repair PSD
    % Inputs:
    %   X: Input covariance matrix
    %   target_max: (Optional) Target maximum value for scaling (default 1000)
    %   verbose: (Optional) Boolean to toggle debug prints (default false)
    
    if nargin < 3, verbose = false; end
    if nargin < 2, target_max = 1000; end
    
    % Ensure Symmetry first
    X = 0.5 * (X + X');
    
    % --- 1. ANALYZE & CLIP RAW OUTLIERS ---
    X_abs = abs(X(:));
    
    % Use 99th percentile for robust max
    raw_robust_max = prctile(X_abs, 99.0); 
    raw_absolute_max = max(X_abs);
    
    % Safety for sparse matrices
    if raw_robust_max == 0, raw_robust_max = raw_absolute_max; end
    if raw_robust_max == 0, raw_robust_max = 1; end 
    
    % Define Threshold (3x the 99th percentile)
    clip_threshold_raw = raw_robust_max * 3.0;
    
    if verbose
        fprintf('  [Outlier Check] Abs Max: %.1f | Threshold: %.1f\n', ...
                raw_absolute_max, clip_threshold_raw);
    end
            
    num_clipped = sum(X_abs > clip_threshold_raw);
    
    if num_clipped > 0
        if verbose
            fprintf('  ⚠️ Clipping %d values (Outliers > 3x P99)\n', num_clipped);
        end
        
        % Hard Clip
        X(X > clip_threshold_raw) = clip_threshold_raw;
        X(X < -clip_threshold_raw) = -clip_threshold_raw;
        
        % Re-calculate robust max after clipping for scaling
        raw_robust_max = prctile(abs(X(:)), 99.0);
    end
    
    % --- 2. SCALE ---
    scale_factor = target_max / raw_robust_max;
    X_scaled = X * scale_factor;
    
    % --- 3. FINAL PSD PROJECTION ---
    X_scaled = 0.5 * (X_scaled + X_scaled'); % Ensure symmetry again
    
    [V, D] = eig(X_scaled);
    d_vals = diag(D);
    
    min_eig = min(d_vals);
    if min_eig < -1e-6 
        if verbose
            fprintf('  ⚠️ Matrix became indefinite after clipping (Min Eig: %.2e). Projecting to PSD...\n', min_eig);
        end
        
        % clip negative eigenvalues
        d_vals(d_vals < 0) = 0;
        
        % Reconstruct
        X_scaled = V * diag(d_vals) * V';
        X_scaled = 0.5 * (X_scaled + X_scaled'); % Final symmetry enforcement
    elseif min_eig < 0
        % Tiny numerical noise fix
        X_scaled = X_scaled + (-min_eig + 1e-9) * eye(size(X_scaled));
    end
    
    % --- 4. FINAL REPORT ---
    if verbose
        scaled_abs = abs(X_scaled(:));
        fprintf('  [Final Stats] P99 < %.1f | Absolute Max = %.1f | Min Eig = %.2e\n', ...
            prctile(scaled_abs, 99.0), max(scaled_abs), min(eig(X_scaled)));
    end
    
    X_clean = X_scaled;
end