% =========================================================================
% DYNAMIC EFM MERGER & EXACT DUPLICATE REMOVAL
% =========================================================================

% 1. Define your undefined number of files here (just keep adding to the list)
file_list = {
    '/Users/benedikthaarscheidt/COVlux/data/models/E_coli/efms/efms_matrix_iML1515_gapfill.mat',
    '/Users/benedikthaarscheidt/COVlux/data/models/E_coli/efms/efms_matrix_iML1515_gapfill2.mat',
    '/Users/benedikthaarscheidt/COVlux/data/models/E_coli/efms/efms_matrix_iML1515_2gapfills.mat',
    '/Users/benedikthaarscheidt/COVlux/data/models/E_coli/efms/efms_matrix_iML1515_forBiomass.mat',
    '/Users/benedikthaarscheidt/COVlux/data/models/E_coli/efms/efms_matrix_iML1515_forBiomass2.mat',
    '/Users/benedikthaarscheidt/COVlux/data/models/E_coli/efms/efms_matrix_iML1515_denovo.mat',
    '/Users/benedikthaarscheidt/COVlux/data/models/E_coli/efms/efms_matrix_iML1515_lastbit.mat',
    '/Users/benedikthaarscheidt/COVlux/data/models/E_coli/efms/efms_matrix_iML1515_repair.mat',
    
};

% Initialize empty containers for the loop
EFM_matrix_all  = [];
EFM_support_all = [];
EFM_anchor_all  = [];
EFM_supps_all   = {};
varNames_all    = {};

fprintf('\n--- Merging %d EFM Files ---\n', length(file_list));

% 2. Loop over every file and concatenate
for i = 1:length(file_list)
    fprintf('Loading file %d: %s\n', i, file_list{i});
    D = load(file_list{i});
    
    % Concatenate matrices column-wise
    EFM_matrix_all  = [EFM_matrix_all,  D.EFM_matrix];
    EFM_support_all = [EFM_support_all, D.EFM_support];
    
    % Concatenate metadata (forced to row vectors)
    EFM_anchor_all  = [EFM_anchor_all,  D.EFM_anchor(:)'];
    EFM_supps_all   = [EFM_supps_all,   D.EFM_supps(:)'];
    varNames_all    = [varNames_all,    D.varNames(:)'];
    
    % Grab static model variables from the FIRST file only
    if i == 1
        rowNames = D.rowNames;
        rxnNames = D.rxnNames;
        S        = D.S;
        model_ir = D.model_ir;
    end
end

total_raw_efms = size(EFM_matrix_all, 2);
fprintf('\nTotal EFMs before duplicate removal: %d\n', total_raw_efms);

% =========================================================================
% 3. THE MAGIC 'UNIQUE' ONE-LINER
% =========================================================================
fprintf('Running fast unique duplicate removal...\n');

% MATLAB's unique() works on ROWS. Because EFMs are COLUMNS, we transpose (').
% We also round to 8 decimals to prevent 1e-16 floating-point noise from 
% tricking MATLAB into thinking two identical EFMs are different.
[~, unique_idx] = unique(round(EFM_matrix_all', 8), 'rows', 'stable');

% 4. Slice EVERYTHING using the safe, unique indices
EFM_matrix  = EFM_matrix_all(:, unique_idx);
EFM_support = EFM_support_all(:, unique_idx);
EFM_anchor  = EFM_anchor_all(unique_idx);
EFM_supps   = EFM_supps_all(unique_idx);

% 5. Fix varNames (column headers) to ensure they are strictly unique
% (If file 1 and 2 both had an "EFM_1_anchor_x", this fixes the name clash)
varNames_raw = varNames_all(unique_idx);
varNames     = matlab.lang.makeUniqueStrings(varNames_raw);

total_clean_efms = size(EFM_matrix, 2);
fprintf('Removed %d exact duplicates.\n', total_raw_efms - total_clean_efms);
fprintf('Final Unique EFMs: %d\n', total_clean_efms);

% =========================================================================
% 6. SAVE MASTER FILE
% =========================================================================
targetFile = '/Users/benedikthaarscheidt/COVlux/data/models/E_coli/efms/efms_matrix_iML1515_MASTER_CLEAN.mat';
fprintf('Saving master basis to %s...\n', targetFile);

save(targetFile, 'EFM_matrix', 'EFM_support', 'EFM_anchor', 'EFM_supps', ...
     'varNames', 'rowNames', 'rxnNames', 'S', 'model_ir', '-v7.3');

fprintf('SUCCESS!  file is ready to use.\n\n');