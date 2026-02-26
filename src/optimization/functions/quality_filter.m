function [clean_mat, mask] = quality_filter(mat, snr, noise)
    residuals = mat(1, :);
    fluxes    = mat(2:end, :);
    n = size(mat, 2);
    mask = true(1, n);
    
    % --- 1. SNR / GHOST FILTER ---
    for k = 1:n
        res = residuals(k);
        v = abs(fluxes(:, k));
        
        active = v(v > noise);
        if isempty(active), min_flux = 0; else, min_flux = min(active); end
        
        if res == 0, ratio = Inf;
        elseif min_flux == 0, ratio = 0;
        else, ratio = min_flux / res;
        end
        
        if ratio < snr
            mask(k) = false;
        end
    end
    fprintf('Removed %d "Ghost" EFMs.\n', sum(~mask));
    
    % --- 2. JACCARD REDUNDANCY FILTER ---
    jaccard_threshold = 0.9;
    valid_idx = find(mask); % EFMs that survived the Ghost filter
    n_valid = length(valid_idx);
    
    if n_valid > 0
        % Extract valid residuals and build boolean support matrix
        valid_residuals = residuals(valid_idx);
        valid_fluxes = fluxes(:, valid_idx);
        supports = abs(valid_fluxes) > noise; 
        
        % SORTING: Prioritize the highest quality EFMs (Lowest Residual Error)
        [~, sort_order] = sort(valid_residuals, 'ascend');
        
        accepted_local_idx = []; 
        
        for i = 1:n_valid
            candidate_idx = sort_order(i);
            candidate_support = supports(:, candidate_idx);
            
            is_redundant = false;
            
            % Vectorized Jaccard check against all currently accepted EFMs
            if ~isempty(accepted_local_idx)
                accepted_supports_mat = supports(:, accepted_local_idx);
                
                % Fast Jaccard Math: Intersection / Union
                intersections = candidate_support' * accepted_supports_mat;
                unions = sum(candidate_support) + sum(accepted_supports_mat, 1) - intersections;
                jaccard_sims = intersections ./ unions;
                
                % If it overlaps > 60% with ANY accepted EFM, flag it
                if any(jaccard_sims > jaccard_threshold)
                    is_redundant = true;
                end
            end
            
            if ~is_redundant
                accepted_local_idx = [accepted_local_idx, candidate_idx];
            else
                % Update global mask to drop this redundant EFM
                global_idx = valid_idx(candidate_idx);
                mask(global_idx) = false;
            end
        end
        fprintf('Removed %d Redundant EFMs (Jaccard > %.1f).\n', n_valid - length(accepted_local_idx), jaccard_threshold);
    end
    
    clean_mat = mat(:, mask);
end
