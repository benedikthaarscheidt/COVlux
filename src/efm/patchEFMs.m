% === Configuration ===
% List all three files here

efmMatPath  = '/Users/benedikthaarscheidt/M.Sc./master_thesis_second_moment/efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5_higher_cover.mat';
efmMatPath2 = '/Users/benedikthaarscheidt/M.Sc./master_thesis_second_moment/efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5.mat';
efmMatPath3 = '/Users/benedikthaarscheidt/M.Sc./master_thesis_second_moment/efms_NEW_iML1515_pruned_permissive_biomassrm_loopless_v5.mat';
filesToPatch = {efmMatPath, efmMatPath2, efmMatPath3}; 

targetName      = 'BIOMASS_Ec_iML1515_WT_75p37M';
targetRxnIdx    = 2219;       % The index in the name list
targetMatrixRow = 2220;       % The index in the matrix (1 + 2219)

% === Master Loop ===
for k = 1:length(filesToPatch)
    fpath = filesToPatch{k};
    fprintf('\n------------------------------------------------\n');
    fprintf('Processing File %d: %s\n', k, fpath);
    
    % 1. Load Data
    data = load(fpath);
    
    % 2. Auto-Detect Variable Names
    % Find the Flux Matrix
    if isfield(data, 'EFM_matrix')
        matVar = 'EFM_matrix';
    elseif isfield(data, 'new_EFM_matrix')
        matVar = 'new_EFM_matrix';
    else
        warning('Skipping %s: No EFM matrix found.', fpath); continue;
    end
    M = data.(matVar);
    
    % Find the Reaction Names
    if isfield(data, 'rxnNames')
        nameVar = 'rxnNames';
    elseif isfield(data, 'rxnE')
        nameVar = 'rxnE';
    else
        warning('Skipping %s: No reaction names found.', fpath); continue;
    end
    names = string(data.(nameVar)(:));
    
    % 3. Check if Patch is Needed
    if any(strcmp(names, targetName))
        fprintf('  [OK] Biomass row already exists here. No changes.\n');
        continue;
    end
    
    % 4. Perform Surgery (Insert Zero Row)
    fprintf('  [PATCHING] Inserting empty Biomass row at #%d...\n', targetMatrixRow);
    
    numCols = size(M, 2);
    zeroRow = zeros(1, numCols);
    
    % Splice Matrix
    M_new = [M(1:targetMatrixRow-1, :); ...
             zeroRow; ...
             M(targetMatrixRow:end, :)];
         
    % Splice Names
    names_new = [names(1:targetRxnIdx-1); ...
                 string(targetName); ...
                 names(targetRxnIdx:end)];
             
    % 5. Save Permanently
    data.(matVar)  = M_new;
    data.(nameVar) = cellstr(names_new); % Ensure it stays as cell array/string
    
    save(fpath, '-struct', 'data');
    fprintf('  [SUCCESS] Saved updated file.\n');
end

fprintf('\n------------------------------------------------\n');
fprintf('All 3 files processed.\n');