% 1. Load both files into separate structures
fprintf('Loading Base Basis (lastbit)...\n');
D1 = load('/Users/benedikthaarscheidt/COVlux/efms_matrix_iML1515_lastbit.mat');

fprintf('Loading Gapfill Basis...\n');
D2 = load('/Users/benedikthaarscheidt/COVlux/efms_matrix_iML1515_gapfill2.mat');

fprintf('Merging datasets...\n');

% 2. Concatenate the massive flux and binary support matrices column-wise
EFM_matrix  = [D1.EFM_matrix, D2.EFM_matrix];
EFM_support = [D1.EFM_support, D2.EFM_support];

% 3. Concatenate the metadata (forcing them to row vectors to be safe)
EFM_anchor  = [D1.EFM_anchor(:)', D2.EFM_anchor(:)']; 
EFM_supps   = [D1.EFM_supps(:)',  D2.EFM_supps(:)'];
varNames    = [D1.varNames(:)',   D2.varNames(:)'];

% 4. Grab the static model variables (they are identical in both files)
rowNames = D1.rowNames;
rxnNames = D1.rxnNames; % This guarantees COVlux line 99 works!
S        = D1.S;
model_ir = D1.model_ir;

% 5. Save the combined master file
targetFile = '/Users/benedikthaarscheidt/COVlux/efms_matrix_iML1515_full_merged.mat';
fprintf('Saving to %s...\n', targetFile);
save(targetFile, 'EFM_matrix', 'EFM_support', 'EFM_anchor', 'EFM_supps', ...
     'varNames', 'rowNames', 'rxnNames', 'S', 'model_ir', '-v7.3');

fprintf('SUCCESS! Merged Basis now has %d EFMs.\n', size(EFM_matrix, 2));