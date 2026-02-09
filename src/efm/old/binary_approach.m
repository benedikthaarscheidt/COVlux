%% binaryEFMs_full.m
% Comprehensive EFM enumeration for irreversible models
clearvars; clc;

%% 1) Load model
fprintf('1) Loading irreversible model...\n');
data = load('e_coli_core_splitandpruned.mat','pruned_ir');
model_ir = data.pruned_ir;
assert(all(model_ir.lb >= 0), 'Model contains reversible reactions!');

%% 2) Enumerate EFMs
fprintf('2) Enumerating EFMs...\n');
[EMbin, EMcoef] = enumerateEFMs_irreversible(model_ir);
fprintf('   Found %d EFMs\n', size(EMbin,2));

%% 3) Validate coverage
covered_rxns = any(EMbin, 2);
uncovered = find(~covered_rxns);
if ~isempty(uncovered)
    fprintf('   Warning: %d reactions not in any EFM:\n', length(uncovered));
    disp(model_ir.rxns(uncovered));
else
    fprintf('   All reactions participate in at least one EFM\n');
end

%% 4) Save results
save('EFMs_binary_full.mat', 'EMbin', 'EMcoef');
fprintf('Done.\n');

%% Core EFM Enumeration Function
function [EFM_binary, EFM_coefficients] = enumerateEFMs_irreversible(model)
    % Robust EFM enumeration with improved numerical handling
    N = full(model.S);  
    [m, q] = size(N);
    
    %% 1. Compute Null Space with Conservative Tolerance
    fprintf('Computing null space... ');
    [U, S, V] = svd(N);
    singular_values = diag(S);
    
    % More conservative tolerance for biological systems
    tol = max(size(N)) * eps(singular_values(1)) * 100; % 100x more permissive
    rankN = sum(singular_values > tol);
    K = V(:, rankN+1:end);
    nullity = size(K,2);
    fprintf('Rank %d, Nullity %d\n', rankN, nullity);
    
    %% 2. Initialize Rays with All Possible Modes
    R1 = (abs(K) > tol);  % Binary mask
    R2 = K;               % Real coefficients
    num_rays = size(R2, 2);
    fprintf('Initialized %d candidate modes\n', num_rays);
    
    %% 3. Enhanced Double-Description Method
    fprintf('Generating EFMs...\n');
    
    for p = 1:q
        fprintf('Processing row %d/%d...\n', p, q);
        
        % Split by current row's sign with relaxed tolerance
        current_row = R2(p,:);
        jneg = find(current_row < -tol/10);  % More permissive
        jpos = find(current_row > tol/10);
        
        % Generate new rays
        new_R1 = false(q, length(jneg)*length(jpos));
        new_R2 = zeros(q, length(jneg)*length(jpos));
        new_count = 0;
        
        for k = jneg
            for l = jpos
                new_support = R1(:,k) | R1(:,l);
                
                % Relaxed minimality check
                is_minimal = true;
                for r = 1:num_rays
                    if r ~= k && r ~= l && all(new_support >= R1(:,r))
                        is_minimal = false;
                        break;
                    end
                end
                
                if is_minimal
                    new_count = new_count + 1;
                    new_R1(:,new_count) = new_support;
                    new_R2(:,new_count) = current_row(l)*R2(:,k) - current_row(k)*R2(:,l);
                end
            end
        end
        
        % Merge new rays
        if new_count > 0
            R1 = [R1(:,setdiff(1:num_rays, jneg)), new_R1(:,1:new_count)];
            R2 = [R2(:,setdiff(1:num_rays, jneg)), new_R2(:,1:new_count)];
            num_rays = size(R2, 2);
            fprintf('  New EFMs: %d (Total: %d)\n', new_count, num_rays);
        end
    end
    
    %% 4. Post-Processing with Biological Validation
    fprintf('Post-processing %d candidate EFMs...\n', num_rays);
    
    % Remove duplicates
    [~, unique_idx] = unique(round(double(R1'),5), 'rows', 'stable'); % Relaxed rounding
    R1 = R1(:,unique_idx);
    R2 = R2(:,unique_idx);
    
    % Biological relevance filter
    valid_efms = true(1, size(R2,2));
    for e = 1:size(R2,2)
        % Check for trivial modes (e.g., single reactions)
        if nnz(R1(:,e)) <= 2
            valid_efms(e) = false;
            continue;
        end
        
        % Check mass balance
        if any(abs(N*R2(:,e)) > 1e-6)
            valid_efms(e) = false;
        end
    end
    
    R1 = R1(:,valid_efms);
    R2 = R2(:,valid_efms);
    
    %% 5. Final Validation and Output
    fprintf('Found %d biologically valid EFMs\n', size(R1,2));
    
    % Reaction participation analysis
    covered_rxns = sum(R1, 2) > 0;
    uncovered = find(~covered_rxns);
    
    if ~isempty(uncovered)
        fprintf('%d reactions not in any EFM:\n', length(uncovered));
        % Display first 20 for brevity
        disp(model.rxns(uncovered(1:min(20,end)))); 
    end
    
    EFM_binary = R1;
    EFM_coefficients = R2;
end