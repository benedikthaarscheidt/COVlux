function [A_opt, selected_idx, admm_metrics, log_Iter, log_K, log_GlobalID] = solve_A_L21_ADMM(E, X, lambda_fraction, max_iter, tol, global_idx, start_K)
    if nargin < 4, max_iter = 1000; end
    if nargin < 5, tol = 1e-4; end
    n_efms = size(E, 2);
    
    % Tracking array: Will store the LAST iteration an EFM was active
    last_active_iter = zeros(n_efms, 1);
    
    % --- 1. PRECOMPUTE CONSTANTS ---
    H = E' * E;
    C_data = E' * X * E;
    
    [V, D] = eig((H + H')/2); 
    d = diag(D);
    
    rho = mean(abs(d));
    if rho < 1e-4 || isnan(rho), rho = 1.0; end 
    Denom = d * d' + rho; 
    
    % --- 2. SCALING ---
    max_signal = max(sqrt(sum(C_data, 2)));
    lambda = lambda_fraction * max_signal; 
    tau_threshold = lambda / rho; 
    
    fprintf('  Starting ADMM L2,1 (Lambda Fraction: %.3f -> Tau: %.4f)...\n', lambda_fraction, tau_threshold);
    
    % --- 3. INITIALIZE VARIABLES ---
    A = zeros(n_efms, n_efms);
    Z = zeros(n_efms, n_efms);
    U = zeros(n_efms, n_efms); 
    
    % --- 4. ADMM LOOP ---
    for k = 1:max_iter
        A_old = A;
        Z_old = Z;
        
        % STEP 1: A-Update
        M = C_data + rho * (Z - U);
        M = (M + M') / 2; 
        
        M_tilde = V' * M * V;
        A_tilde = M_tilde ./ Denom;
        A = V * A_tilde * V';
        A = (A + A') / 2; 
        
        % STEP 2: Z-Update
        W = A + U;
        Z = zeros(n_efms, n_efms);
        row_norms = sqrt(sum(W.^2, 2));
        
        % Update true active tracking
        active_now = row_norms > tau_threshold;
        last_active_iter(active_now) = k; 
        
        %for i = 1:n_efms
        %    if active_now(i)
        %        Z(i, :) = W(i, :) * (1 - tau_threshold / row_norms(i));
        %    else
        %        Z(i, :) = 0; % Nuke to zero
        %    end
        %end
        
        c = zeros(n_efms, 1);
        c(active_now) = 1 - (tau_threshold ./ row_norms(active_now));
        
        
        Z = W .* (c * c');
        
        
        U = U + A - Z;
        
        % CONVERGENCE CHECK
        r_norm = norm(A - Z, 'fro'); 
        s_norm = norm(-rho * (Z - Z_old), 'fro'); 
        
        if isnan(r_norm) || isinf(r_norm)
            error('ADMM destabilized! Check data scaling.');
        end
        
        if mod(k, 50) == 0 || k == 1
            current_survivors = sum(active_now);
            fprintf('    Iter %4d | Primal: %10.4f | Dual: %10.4f | EFMs Kept: %d / %d\n', ...
                k, r_norm, s_norm, current_survivors, n_efms);
        end
        
        if r_norm < tol && s_norm < tol
            current_survivors = sum(active_now);
            fprintf('  -> ADMM converged at iter %d. Final EFMs Kept: %d\n', k, current_survivors);
            break;
        end
    end
    
    if k == max_iter
        fprintf('  -> ADMM reached max iterations (%d)\n', max_iter);
    end
    
    k_final = k;
    
    % --- 5. EXTRACT SURVIVORS AND BUILD LOG ---
    % EFMs that survived all the way to k_final
    survivors = (last_active_iter == k_final);
    selected_idx = find(survivors);
    dropped_local_idx = find(~survivors);
    
    % Sort the dropped EFMs by when they died (ascending)
    [death_iters, sort_order] = sort(last_active_iter(dropped_local_idx), 'ascend');
    dropped_sorted_local = dropped_local_idx(sort_order);
    
    % 5a. Build logs for the DROPPED EFMs
    log_Iter_dropped = death_iters + 1; % "Died" the iteration after it was last active
    log_GlobalID_dropped = global_idx(dropped_sorted_local);
    log_GlobalID_dropped = log_GlobalID_dropped(:);
    log_Iter_dropped = log_Iter_dropped(:);
    
    log_K_dropped = zeros(length(log_GlobalID_dropped), 1);
    curr_K = start_K;
    for i = 1:length(log_GlobalID_dropped)
        curr_K = curr_K - 1;
        log_K_dropped(i) = curr_K;
    end
    
    % 5b. Build dummy logs for the SURVIVORS (So the biological script sees them)
    survivor_global_ids = global_idx(selected_idx);
    log_Iter_survivors = repmat(k_final + 1, length(survivor_global_ids), 1);
    log_GlobalID_survivors = survivor_global_ids(:);
    log_K_survivors = zeros(length(survivor_global_ids), 1);
    
    % Count K down to 0 for the survivors
    for i = 1:length(log_GlobalID_survivors)
        curr_K = curr_K - 1;
        log_K_survivors(i) = curr_K;
    end
    
    % Combine them
    log_Iter = [log_Iter_dropped; log_Iter_survivors];
    log_K = [log_K_dropped; log_K_survivors];
    log_GlobalID = [log_GlobalID_dropped; log_GlobalID_survivors];
    
    A_opt = A;
    admm_metrics.iterations = k;
    admm_metrics.primal_residual = r_norm;
    admm_metrics.efms_kept = length(selected_idx);
end