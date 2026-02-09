function [num_unique_lost, lost_names, final_efm_idx, revived_names] = ...
    compute_unique_reaction_loss(A_opt, E_calc, mu_vector, E_target, rxn_names_target,verbose)
% COMPUTE_UNIQUE_REACTION_LOSS
% Updated for Second Moment Pipeline: "Rescue" logic is removed.
%
% Rationale: 
%   Since Target X = Covariance + Lambda * (Mean*Mean'), the solver 
%   is mathematically forced to select high-mean EFMs to fit the diagonal.
%   Therefore, 'active_idx' (Lasso) already includes Housekeeping pathways.
%
% Inputs:
%   A_opt:      Optimized covariance matrix (The Truth).
%   E_calc:     (Unused) Kept for compatibility.
%   mu_vector:  (Unused) Kept for compatibility.
%   E_target:   Matrix to calculate loss on (Full or Reduced).
%   rxn_names_target: Reaction names for E_target.

    [n_target_rxns, ~] = size(E_target);
    
    % --- 1. IDENTIFY ACTIVE EFMS ---
    % Trust the solver. Any EFM with non-zero variance contributes 
    % to explaining either the Fluctuations or the Mean Flux.
    variances = diag(A_opt);
    final_efm_idx = find(variances > 1e-9);
    
    % --- 2. CALCULATE LOSS ON TARGET MATRIX ---
    if ~isempty(final_efm_idx)
        % Map selected EFMs to reactions in the TARGET space
        % Check which reactions have flux > 1e-6 in at least one selected EFM
        reactions_kept_mask = any(abs(E_target(:, final_efm_idx)) > 1e-6, 2);
    else
        reactions_kept_mask = false(n_target_rxns, 1);
    end
    
    num_kept = sum(reactions_kept_mask);
    lost_names = rxn_names_target(~reactions_kept_mask);
    num_unique_lost = length(lost_names);
    
    % --- 3. REVIVAL IS OBSOLETE ---
    % We return empty to maintain function signature compatibility
    revived_names = []; 
    
    if verbose
    % --- 4. LOGGING ---
        fprintf('  [Loss Analysis] Total EFMs Selected: %d\n', length(final_efm_idx));
        fprintf('  [Loss Analysis] Reactions Kept: %d (%.1f%%) | Lost: %d\n', ...
            num_kept, (num_kept/n_target_rxns)*100, num_unique_lost);
    end 
end
