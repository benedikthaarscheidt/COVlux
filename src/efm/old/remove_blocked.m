%% remove_blocked_and_split_biomass.m

cobraModelFile = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/iML1515.mat';

%% 1) Load original reversible model
orig = readCbModel(cobraModelFile);
fprintf('1) Original reaction count (reversible iML1515): %d\n', numel(orig.rxns));
tol = 9e-6;
%% 2) Convert to irreversible
model_ir = convertToIrreversible(orig);
fprintf('2) Converted to irreversible: %d reactions.\n\n', numel(model_ir.rxns));


biomass_idx = find(contains(model_ir.rxns, "BIOMASS", "IgnoreCase", true), 1);
if isempty(biomass_idx)
    error('No “BIOMASS” reaction found in the pruned irreversible model.');
end

% Extract its stoichiometric column
bio_col = model_ir.S(:, biomass_idx);
consumed_idx     = find(bio_col <  0);   % rows where biomass pulled in these metabolites
produced_idx     = find(bio_col >  0);   % rows where biomass produced these metabolites
consumed_coeffs  = bio_col(consumed_idx); 



%% 5) Add PRO_<met> (source) and DM_<met> (sink) reactions for each biomass metabolite

% Record current reaction count so we know which ones are newly added
[~, nOrig] = size(model_ir.S);

% 5.A) For each consumed metabolite, add a sink DM_<met>
for p = 1:numel(consumed_idx)
    metRow = consumed_idx(p);
    metID  = model_ir.mets{metRow};  
    coeff  = consumed_coeffs(p);      

    % Reaction name “DM_<safeMetID>”
    safeMetID = regexprep(metID, '[^A-Za-z0-9_]', '_');
    rxnName   = sprintf('DM_%s', safeMetID);

    % Build new column: 0 everywhere except –1 at metRow
    m_old   = size(model_ir.S,1);
    newCol  = zeros(m_old,1);
    newCol(metRow) = -1;

    model_ir.S       = [model_ir.S, newCol];
    model_ir.rxns{end+1} = rxnName;
    model_ir.lb(end+1)   = 0;      % irreversible sink
    model_ir.ub(end+1)   = Inf;   % large UB

    fprintf('   Added %s:  %g·%s → ∅\n', rxnName, coeff, metID);
end

fprintf('\n   Removing original biomass reaction “%s” (column %d) …\n', ...
        model_ir.rxns{biomass_idx}, biomass_idx);
model_ir.S(:,      biomass_idx) = [];
model_ir.rxns(    biomass_idx) = [];
model_ir.lb(      biomass_idx) = [];
model_ir.ub(      biomass_idx) = [];
fprintf('   Biomass reaction removed.\n\n');

%% 3) Iteratively remove blocked reactions from the irreversible model

iteration = 0;

while true
    iteration    = iteration + 1;
    removedCount = 0;
    fprintf('3.%d) Scanning for blocked reactions …\n', iteration);
    
    % Start at the first reaction and work your way through
    j = 1;
    while j <= numel(model_ir.rxns)
        % 1) Build a temp model that maximizes v_j
        tmp      = model_ir;
        nRxns    = numel(tmp.rxns);
        tmp.c    = zeros(nRxns,1);
        tmp.c(j) = 1;
        
        sol_max = optimizeCbModel(tmp, 'max');
        
        % 2) If it can't carry flux ≳ tol, it's blocked
        if isempty(sol_max.x) || sol_max.x(j) < tol
            fprintf('   → Removing blocked reaction %s (index %d)\n', ...
                    model_ir.rxns{j}, j);
            model_ir = removeRxns(model_ir, model_ir.rxns{j}, 0, 1);
            removedCount = removedCount + 1;
            % ! DO NOT increment j: the next reaction has just shifted into position j
        else
            % Reaction j is OK: move on
            j = j + 1;
        end
    end
    
    % 3) If nothing was removed this pass, we’re done
    if removedCount == 0
        fprintf('3.%d) No more blocked reactions. Stopping prune.\n\n', iteration);
        break;
    else
        fprintf('   ► Removed %d reactions this pass; %d remain in model.\n\n', ...
                removedCount, numel(model_ir.rxns));
    end
end

pruned_ir = model_ir;  % this is our fully‐pruned irreversible model
fprintf('Final pruned irreversible model has %d reactions.\n\n', numel(pruned_ir.rxns));

disp(iteration)

%%         Uncomment below if you also need PRO_ (otherwise downstream may remain starved).
%for p = 1:numel(produced_idx)
%    metRow = produced_idx(p);
%    metID  = pruned_ir.mets{metRow};
%    coeff  = produced_coeffs(p);
%
%    safeMetID = regexprep(metID, '[^A-Za-z0-9_]', '_');
%    rxnName   = sprintf('PRO_%s', safeMetID);
%
%    m_old   = size(pruned_ir.S,1);
%    newCol  = zeros(m_old,1);
%    newCol(metRow) =  coeff;
%
%    pruned_ir.S       = [pruned_ir.S, newCol];
%    pruned_ir.rxns{end+1} = rxnName;
%    pruned_ir.lb(end+1)   = 0;      % irreversible source
%    pruned_ir.ub(end+1)   = 1000;   % large UB
%
%    fprintf('   Added %s:  ∅ → %g·%s\n', rxnName, coeff, metID);
%end
%
% 5.C) Remove the original biomass reaction itself

r = 1015;
eps_test = 9e-6;        % same as your eps_flux

% 1) Copy pruned model and force v(r) ≥ eps_test
tmp = model_ir;         
tmp.lb(r) = eps_test;   % force r active at at least eps_test
tmp.ub(r) = model_ir.ub(r);
tmp.c   = zeros(length(tmp.rxns),1);

% 2) Solve a zero‐objective FBA (just check feasibility)
sol = optimizeCbModel(tmp, 'max');

if isempty(sol.x) || sol.stat ~= 1
    fprintf('\n→ [CHECK 1] Reaction %d (%s) is already blocked in pruned_ir: no steady‐state flux possible.\n\n', ...
            r, model_ir.rxns{r});
    return
else
    fprintf('\n→ [CHECK 1] OK: There _is_ a v-vector with v(%d) = %.2e.\n\n', ...
            r, sol.x(r));
end


%% 7) Save the final pruned + split‐biomass irreversible model
outFile = '/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/pruned_iML1515_irrev_split_biomassrm.mat';
save(outFile, 'pruned_ir');
fprintf('7) Final pruned+split irreversible model saved to:\n   %s\n', outFile);
fprintf('   Total reactions now: %d\n', numel(pruned_ir.rxns));
