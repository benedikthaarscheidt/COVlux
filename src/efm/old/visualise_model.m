clc; clear;
data     = load('/Users/benedikthaarscheidt/M.Sc./master_thesis/Models/generic_models/E_coli/e_coli_core_splitandpruned.mat','pruned_ir');
model_ir = data.pruned_ir;
S        = sparse(model_ir.S);
rxnNames = cellstr(model_ir.rxns(:));
[m,n]    = size(S);
eps_flux    = 1e-5;
tol_balance = 1e-6;
Mdefault    = 1e3;
if isfield(model_ir,'ub'), ubv = model_ir.ub(:); else, ubv = Mdefault*ones(n,1); end
if isempty(ubv), ubv = Mdefault*ones(n,1); end
ubv(~isfinite(ubv)) = Mdefault;
if isfield(model_ir,'lb'), lbv = max(0,model_ir.lb(:)); else, lbv = zeros(n,1); end
lb  = [lbv; zeros(n,1)];
ub  = [ubv;  ones(n,1)];
vtype = [ repmat('C',n,1); repmat('B',n,1)];
Aeq = [S, sparse(m,n)];
beq = zeros(m,1);
A_link = [speye(n), -spdiags(ubv,0,n,n); -speye(n), eps_flux*speye(n)];
b_link = zeros(2*n,1);
onlyProduces = all(S>=0,1)' & any(S>0,1)';
onlyConsumes = all(S<=0,1)' & any(S<0,1)';
A_imp = []; b_imp = [];
if any(onlyProduces)
    r = sparse(1, 2*n); r(1, n + find(onlyProduces)) = -1; A_imp = r; b_imp = -1;
end
A_exp = []; b_exp = [];
if any(onlyConsumes)
    r = sparse(1, 2*n); r(1, n + find(onlyConsumes)) = -1; A_exp = r; b_exp = -1;
end
Aineq = [A_link; A_imp; A_exp];
bineq = [b_link; b_imp; b_exp];
mdlBase.A          = [Aeq; Aineq];
mdlBase.rhs        = [beq; bineq];
mdlBase.sense      = [repmat('=',size(Aeq,1),1); repmat('<',size(Aineq,1),1)];
mdlBase.lb         = lb;
mdlBase.ub         = ub;
mdlBase.vtype      = vtype;
mdlBase.modelsense = 'min';
mdlBase.obj        = [zeros(n,1); ones(n,1)];
params = struct('OutputFlag',0,'FeasibilityTol',1e-9,'IntFeasTol',1e-9,'NumericFocus',3,'TimeLimit',30);
badPartner = zeros(n,1);
for i = 1:n
    name_i = rxnNames{i};
    if endsWith(name_i,'_f')
        baseName = name_i(1:end-2);
        j = find(strcmp(rxnNames, [baseName '_b']), 1);
        if ~isempty(j)
            if max(abs(S(:,i) + S(:,j))) <= 1e-12
                badPartner(i) = j; badPartner(j) = i;
            end
        end
    end
end
bpairs = find(badPartner>0 & (1:n)' < badPartner);
p      = numel(bpairs);
fprintf('Found %d reversible pairs.\n', p);
results = table('Size',[p 12], ...
    'VariableTypes', {'string','string','logical','logical','logical','logical','logical','string','logical','logical','double','double'}, ...
    'VariableNames', {'rxn_f','rxn_b','f_alone_ok','b_alone_ok','b_needs_twin','f_needs_twin','both_on_ok','class','f_support_is_efm','b_support_is_efm','nnz_f','nnz_b'});
for t = 1:p
    i = bpairs(t); j = badPartner(i);
    if endsWith(rxnNames{j},'_f') && ~endsWith(rxnNames{i},'_f')
        tmp=i; i=j; j=tmp;
    end
    results.rxn_f(t) = string(rxnNames{i});
    results.rxn_b(t) = string(rxnNames{j});
    [feasF, vF] = solve_one(mdlBase, params, n, i, j);
    results.f_alone_ok(t) = feasF;
    if feasF
        suppF = find(vF > max(0,eps_flux - tol_balance));
        results.f_support_is_efm(t) = nullity1_ok(S, suppF);
        results.nnz_f(t) = numel(suppF);
    end
    [feasB, vB] = solve_one(mdlBase, params, n, j, i);
    results.b_alone_ok(t) = feasB;
    if feasB
        suppB = find(vB > max(0,eps_flux - tol_balance));
        results.b_support_is_efm(t) = nullity1_ok(S, suppB);
        results.nnz_b(t) = numel(suppB);
    end
    bNeeds = false; fNeeds = false;
    if ~feasB
        [feasB2, usedTwinB] = solve_allow_twin(mdlBase, params, n, j, i);
        bNeeds = feasB2 && usedTwinB;
    end
    if ~feasF
        [feasF2, usedTwinF] = solve_allow_twin(mdlBase, params, n, i, j);
        fNeeds = feasF2 && usedTwinF;
    end
    MdlBoth = mdlBase;
    MdlBoth.lb(n+i)=1; MdlBoth.ub(n+i)=1;
    MdlBoth.lb(n+j)=1; MdlBoth.ub(n+j)=1;
    solBoth = gurobi(MdlBoth, params);
    bothOK  = isfield(solBoth,'status') && strcmp(solBoth.status,'OPTIMAL');
    results.b_needs_twin(t) = bNeeds;
    results.f_needs_twin(t) = fNeeds;
    results.both_on_ok(t)   = bothOK;
    if results.f_alone_ok(t) && results.b_alone_ok(t)
        results.class(t) = "BOTH_ALONE_OK";
    elseif results.f_alone_ok(t) && ~results.b_alone_ok(t)
        if bNeeds
            results.class(t) = "BACKWARD_NEEDS_TWIN";
        elseif bothOK
            results.class(t) = "BACKWARD_BLOCKED_BUT_BOTH_OK";
        else
            results.class(t) = "BACKWARD_BLOCKED";
        end
    elseif ~results.f_alone_ok(t) && results.b_alone_ok(t)
        if fNeeds
            results.class(t) = "FORWARD_NEEDS_TWIN";
        elseif bothOK
            results.class(t) = "FORWARD_BLOCKED_BUT_BOTH_OK";
        else
            results.class(t) = "FORWARD_BLOCKED";
        end
    else
        if bothOK
            results.class(t) = "BOTH_NEED_EACH_OTHER";
        else
            results.class(t) = "BOTH_BLOCKED";
        end
    end
end
disp(results(:,{'rxn_f','rxn_b','class','f_alone_ok','b_alone_ok','f_needs_twin','b_needs_twin','both_on_ok'}));
fprintf('\nClass counts:\n');
cats = categories(categorical(results.class));
for k = 1:numel(cats)
    c = cats{k};
    fprintf('  %-26s : %d\n', c, sum(results.class==c));
end
function [feas, v] = solve_one(mdlBase, params, n, r_on, r_off)
    Mdl = mdlBase;
    Mdl.lb(n + r_on)  = 1;  Mdl.ub(n + r_on)  = 1;
    Mdl.lb(n + r_off) = 0;  Mdl.ub(n + r_off) = 0;
    sol  = gurobi(Mdl, params);
    feas = isfield(sol,'status') && strcmp(sol.status,'OPTIMAL');
    if feas
        x = sol.x; v = x(1:n);
    else
        v = [];
    end
end
function [feas, usedTwin] = solve_allow_twin(mdlBase, params, n, r_on, twin)
    Mdl = mdlBase;
    Mdl.lb(n + r_on) = 1;  Mdl.ub(n + r_on) = 1;
    sol  = gurobi(Mdl, params);
    feas = isfield(sol,'status') && strcmp(sol.status,'OPTIMAL');
    usedTwin = false;
    if feas
        b = sol.x(n+1:end);
        usedTwin = (b(twin) > 0.5);
    end
end
function isEFM = nullity1_ok(Smat, supp)
    if isempty(supp), isEFM=false; return; end
    A = full(Smat(:,supp));
    s = svd(A,'econ');
    if isempty(s), isEFM=false; return; end
    tol = max(size(A)) * eps(max(s));
    rk  = sum(s > tol);
    isEFM = (numel(supp) - rk == 1);
end
