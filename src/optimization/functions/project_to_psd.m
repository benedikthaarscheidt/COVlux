function X_psd = project_to_psd(X, min_eig)
% Project symmetric matrix to nearest PSD matrix
    if nargin < 2, min_eig = 1e-8; end

    % Ensure symmetry
    X = 0.5 * (X + X');

    % Eigen decomposition
    [V, D] = eig(X);
    eig_vals = diag(D);

    % Clip negative eigenvalues to min_eig
    eig_vals_clean = max(eig_vals, min_eig);

    % Reconstruct
    X_psd = V * diag(eig_vals_clean) * V';

    % Ensure symmetry (numerical)
    X_psd = 0.5 * (X_psd + X_psd');

    % Optional: Verify
    [~, p] = chol(X_psd);
    if p > 0
        % warning('PSD projection failed, adding more regularization');
        X_psd = X_psd + (min_eig * 10) * eye(size(X_psd));
    end
end