function A_final = solve_cvx(E_sel, M)
% SOLVE_FINAL_WEIGHTS_CVX
% Refines the covariance matrix A for a selected subset of EFMs using Convex Optimization.
%
% OBJECTIVE: Minimize || M - E * A * E' ||_F
% CONSTRAINT: A must be Positive Semi-Definite (PSD) and Symmetric.
%
% Inputs:
%   E_sel : The subset of normalized EFMs (n_rxn x k)
%   M     : The target Second Moment / Covariance Matrix (n_rxn x n_rxn)
%
% Output:
%   A_final : The optimized k x k covariance matrix of EFM coefficients.

    [n_rxn, k] = size(E_sel);
    
    nugget = 1e-4;

    fprintf('  > CVX Calibration: Optimizing weights for %d EFMs...\n', k);
    
    
    if ~exist('cvx_begin', 'file')
        warning('CVX is not installed or not in the path. Fallback to Algebraic Projection.');
        A_final = fallback_projection(E_sel, M);
        return;
    end

    try
        cvx_begin quiet
            variable A(k, k) symmetric
            
            % The Objective: Minimize reconstruction error (Frobenius Norm)
            minimize( norm( M - E_sel * A * E_sel', 'fro' ) + nugget * norm(A, 'fro') )
            
            subject to
                % The Constraint: Variance must be non-negative (PSD)
                A == semidefinite(k);
                A >= 0;
        cvx_end
        
        if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
            A_final = full(A);
            fprintf('    [Success] CVX Converged. Optimal Error: %.4e\n', cvx_optval);
        else
            warning('CVX Status: %s. Using fallback projection.', cvx_status);
            A_final = fallback_projection(E_sel, M);
        end
        
    catch ME
        warning('CVX Error: %s. Using fallback projection.', ME.message);
        A_final = fallback_projection(E_sel, M);
    end
end

function A_proj = fallback_projection(E, M)
% Fallback: Algebraic Projection + Eigenvalue Clipping
    lambda = 1e-6;
    k = size(E, 2);
    
    % Ridge Regression Projection (E_dagger * M * E_dagger')
    P = (E' * E + lambda * eye(k)) \ E';
    A_raw = P * M * P';
    
    % Force Symmetry
    A_raw = 0.5 * (A_raw + A_raw');
    
    % Clip Negative Eigenvalues
    [V, D] = eig(A_raw);
    d = diag(D);
    d(d < 1e-9) = 1e-9;
    
    A_proj = V * diag(d) * V';
    fprintf('    [Fallback] Computed Algebraic Projection.\n');
end