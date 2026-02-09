function A_psd = robust_psd_projection(A, min_eig_ratio, verbose)
    % More robust PSD projection
    % verbose: (Optional) boolean to toggle debug prints. Default is false.
    
    if nargin < 3, verbose = false; end
    if nargin < 2
        min_eig_ratio = 1e-12;  % More conservative threshold
    end
    
    % Ensure perfect symmetry with multiple iterations
    for i = 1:3
        A = 0.5 * (A + A');
    end
    
    % Eigen decomposition with better conditioning
    [V, D] = eig(A);
    eigenvals = diag(D);
    
    if verbose
        fprintf('  Eigenvalue range before clipping: %.2e to %.2e\n', min(eigenvals), max(eigenvals));
        fprintf('  Number of Negative eigenvalues: %d\n', sum(eigenvals < 0));
    end
    
    % More aggressive clipping with relative threshold
    max_eig = max(abs(eigenvals));
    threshold = max(min_eig_ratio * max_eig, 1e-12);  % Ensure minimum threshold
    eigenvals_clean = max(eigenvals, threshold);
    
    if verbose
        fprintf('  Eigenvalue range after: %.2e to %.2e\n', min(eigenvals_clean), max(eigenvals_clean));
        fprintf('  Clipping threshold: %.2e\n', threshold);
    end
    
    % Reconstruct with symmetry enforcement
    A_psd = V * diag(eigenvals_clean) * V';
    
    % Multiple symmetry enforcements
    for i = 1:3
        A_psd = 0.5 * (A_psd + A_psd');
    end
    
    % Final small regularization to ensure numerical stability
    [~, p] = chol(A_psd);
    if p > 0
        if verbose
            fprintf('  Adding final regularization (p=%d)\n', p);
        end
        A_psd = A_psd + (threshold * 10) * eye(size(A_psd));
    end
end