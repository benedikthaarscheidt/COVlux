function [clean_mat, mask] = quality_filter(mat, snr, noise)
    residuals = mat(1, :);
    fluxes    = mat(2:end, :);
    n = size(mat, 2);
    mask = true(1, n);
    
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
    clean_mat = mat(:, mask);
end
