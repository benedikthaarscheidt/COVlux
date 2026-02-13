function A_opt = solve_QR_regularized_new(E, X, lambda)
   % ROBUST SOLVER: Truncated SVD Pseudoinverse
    % Prevents explosion by ignoring singular values < tolerance
    
   
    
    % 1. Compute SVD of E
    [U, S, V] = svd(E, 'econ');
    s_vals = diag(S);
    
    % 2. Determine safe tolerance
    
    max_s = max(s_vals);
    tolerance = 1e-6 * max_s; 
    
    % 3. Invert significant singular values only
    keep_idx = s_vals > tolerance;
    s_inv = zeros(size(s_vals));
    s_inv(keep_idx) = 1 ./ s_vals(keep_idx);
    
    %fprintf('    SVD Stats: Max sigma=%.2e, Min kept=%.2e. Discarded %d/%d dims.\n', ...
    %    max_s, min(s_vals(keep_idx)), sum(~keep_idx), length(s_vals));
    
    % 4. Construct Robust Pseudoinverse
   
    E_dagger = V * diag(s_inv) * U';
    
    % 5. Compute A = E_dagger * X * E_dagger'
    A_opt = E_dagger * X * E_dagger';
    
    % 6. Ensure Symmetry
    A_opt = 0.5 * (A_opt + A_opt');

    
    % --- VALIDATION ---
    X_recon = E * A_opt * E';
    

    %err_val = norm(X_recon - X, 'fro') / norm(X, 'fro');
    %fprintf('  Initialization Error: %.6f\n', err_val);
    %fprintf('  Element Range: %.2e to %.2e\n', min(A_opt(:)), max(A_opt(:)));
end
