matFilePath = '/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/EFMs_by_reaction.mat';

% Load only the four variables into a struct
data = load(matFilePath, 'efms', 'covered', 'notCovered', 'finalNotCovered');

% Extract each variable from the struct
efms           = data.efms;           % EFM matrix
covered        = data.covered;        % covered array
notCovered     = data.notCovered;     % notCovered array
finalNotCovered= data.finalNotCovered;% finalNotCovered array

% (Optional) Check their sizes / contents
disp('Size of efms:');        disp(size(efms));
disp('Size of covered:');     disp(size(covered));
disp('Size of notCovered:');  disp(size(notCovered));
disp('Size of finalNotCovered:'); disp(finalNotCovered);