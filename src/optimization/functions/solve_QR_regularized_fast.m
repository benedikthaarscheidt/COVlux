function A_opt = solve_QR_regularized_fast(E, X, lambda)
    [n_rxn, k] = size(E);
    [Q, R] = qr(E, 0); 
    
   
    Gram = R' * R;
    
    reg_term = lambda * eye(k); 
    
    % Standard Ridge Solver
    Gram_Reg = Gram + reg_term;
    
    % Sandwich Projection
    Projected_Variance = E' * X * E;
    
    % Solve: A = inv(Gram_Reg) * Projected_Variance * inv(Gram_Reg)
    A_opt = Gram_Reg \ Projected_Variance / Gram_Reg;
    A_opt = 0.5 * (A_opt + A_opt');
    
    X_recon = E * A_opt * E';
    err_total = norm(X_recon - X, 'fro') / norm(X, 'fro');
    var_exp = 1 - err_total^2; 
end