
fprintf('1) Loading irreversible model...\n');
data = load('e_coli_core_splitandpruned.mat','pruned_ir');
model_ir = data.pruned_ir;
rxnNames=model_ir.rxnNames;
S=model_ir.S;
[m,n]=size(S)
%% 

n = size(S,2);
blocked = false(n,1);

options = optimoptions('linprog','Display','none');

for i = 1:n
    f = zeros(n,1); f(i) = 1;
    [~,~,exitflag] = linprog(-f,[],[],S,model_ir.b,model_ir.lb,model_ir.ub,options);
    if exitflag ~= 1
        blocked(i) = true;
    end
end

blocked_rxns = find(blocked);
disp(['Number of blocked reactions: ', num2str(length(blocked_rxns))]);

if exist('rxnNames','var') && ~isempty(rxnNames)
    disp('Blocked reactions (index and name):');
    for j = 1:length(blocked_rxns)
        fprintf('%4d: %s\n', blocked_rxns(j), rxnNames{blocked_rxns(j)});
    end
else
    disp('Blocked reaction indices:');
    disp(blocked_rxns');
end

%% 


splitPairs = identify_split_pairs(rxnNames);

EFMs = enumerate_efms(S, splitPairs);


function splitPairs = identify_split_pairs(rxnNames)
    % rxnNames: cell array of reaction names
    n = numel(rxnNames);
    splitPairs = [];
    baseNames = cell(n,1);
    fwdIdx = [];
    bwdIdx = [];
    for i = 1:n
        name = rxnNames{i};
        if endsWith(name, '_f')
            base = extractBefore(name, strlength(name)-1);
            baseNames{i} = base;
            fwdIdx = [fwdIdx; i];
        elseif endsWith(name, '_b')
            base = extractBefore(name, strlength(name)-1);
            baseNames{i} = base;
            bwdIdx = [bwdIdx; i];
        end
    end
    for i = 1:numel(fwdIdx)
        base = baseNames{fwdIdx(i)};
        idxB = find(strcmp(baseNames, base) & ismember(1:n, bwdIdx));
        if ~isempty(idxB)
            splitPairs = [splitPairs; fwdIdx(i), idxB(1)];
        end
    end
end

function EFMs = enumerate_efms(S, splitPairs)
    % S: Stoichiometric matrix (m x r)
    % splitPairs: nPairs x 2 array, each row contains [col_pos, col_neg] for a split reversible reaction
    % Returns: EFMs, r x N matrix of EFMs (each column is an EFM)
    
    % Enumerate all extreme rays using double description
    disp('Running Double Description method...');
    EFMs_all = double_description(S);
    disp(['Number of extreme rays found: ', num2str(size(EFMs_all,2))]);

    % Filter out futile cycles: EFMs supported only by both split parts of a reversible reaction
    keep = true(1, size(EFMs_all,2));
    for i = 1:size(EFMs_all,2)
        supp = find(EFMs_all(:,i) > 1e-12); % indices of nonzero reactions
        for j = 1:size(splitPairs,1)
            pair = sort(splitPairs(j,:));   % ensure ascending order
            if isequal(sort(supp), pair)
                keep(i) = false;            % support is only the split pair => futile cycle
                break;
            end
        end
    end
    EFMs = EFMs_all(:, keep);
    disp(['Number of EFMs after futile cycle removal: ', num2str(size(EFMs,2))]);

    % Optional: check that each reaction appears in at least one EFM
    present = any(EFMs > 1e-12, 2);
    missing = find(~present);
    if isempty(missing)
        disp('All reactions appear in at least one EFM.');
    else
        disp('WARNING: Some reactions are not present in any EFM:');
        disp(numel(missing));
    end
end


function R = double_description(A)
    % Input:  A - m x d matrix, defines cone P = {x | Ax=0, x>=0}
    % Output: R - d x n matrix, columns are extreme rays of P

    % 1. Find basis of kernel of A (i.e., N(A)), use row-echelon for numerical stability
    K = null(A, 'r');     % each column is a kernel vector
    R = K;                % Initialize R with the basis of N(A)
    d = size(A, 2);       % dimension

    rho = [];             % Set of indices already considered for non-negativity

    while length(rho) < d
        % 2. Pick next non-negativity constraint to process
        unconsidered = setdiff(1:d, rho);
        j = unconsidered(1);   % You can implement smarter choice if needed

        % 3. Classify rays according to sign of j-th coordinate
        Rj = R(j, :);
        tau_pos = find(Rj > 0);
        tau_zero = find(Rj == 0);
        tau_neg = find(Rj < 0);

        % 4. Find adjacent pairs (this is a simple adjacency test: here, two rays are adjacent if their supports overlap in all but one position, this is simplified for generic position)
        tau_adj = [];
        for i = tau_pos
            for k = tau_neg
                support_i = find(R(:, i));
                support_k = find(R(:, k));
                if sum((R(:, i) ~= 0) & (R(:, k) ~= 0)) >= d-1
                    tau_adj = [tau_adj; i, k];
                end
            end
        end

        % 5. Create new rays
        Rnew = [];
        for idx = 1:size(tau_adj,1)
            i = tau_adj(idx,1);
            k = tau_adj(idx,2);
            p = R(:,i);
            q = R(:,k);
            r_ik = q(j) * p - p(j) * q;
            % Normalize to avoid duplicates (optional)
            if any(r_ik)
                r_ik = r_ik / norm(r_ik);  % normalize
                % Check for duplicates
                if isempty(Rnew) || ~ismembertol(r_ik', Rnew', 1e-8, 'ByRows', true)
                    Rnew = [Rnew, r_ik];
                end
            end
        end

        % 6. Update R: keep those with positive or zero j-th coordinate
        R_pos = R(:, tau_pos);
        R_zero = R(:, tau_zero);
        R = [R_pos, R_zero, Rnew];

        % 7. Update set of considered constraints
        rho = [rho, j];
    end
end