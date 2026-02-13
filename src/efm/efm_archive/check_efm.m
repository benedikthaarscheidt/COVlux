function analyze_efm_set(matfile)
if nargin<1 || isempty(matfile)
    if exist('efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5_higher_cover.mat','file')
        matfile = 'efms_matrix_iML1515_pruned_permissive_biomassrm_loopless_v5_higher_cover.mat';
    else
        [f,p] = uigetfile('*.mat','Select EFM .mat file');
        if isequal(f,0), error('No file selected.'); end
        matfile = fullfile(p,f);
    end
end

outdir = fullfile(pwd,'efm_analysis');
if ~exist(outdir,'dir'), mkdir(outdir); end

load(matfile,'EFM_matrix','EFM_table','rowNames','varNames','EFM_support','efmAnchors','rxnNames','S','model_ir');

eps_flux = 1e-5;
tol_balance = 1e-6;
tol_nz = eps_flux - tol_balance;

n = numel(rxnNames);
if ~exist('EFM_matrix','var') || isempty(EFM_matrix), error('EFM_matrix missing or empty.'); end
numEFMs = size(EFM_matrix,2);
residual = EFM_matrix(1,:).';
V = EFM_matrix(2:end,:);

active = abs(V) > tol_nz;
efm_len = sum(active,1).';
rxn_hits = sum(active,2);
rxn_freq = rxn_hits/numEFMs;

rxn_isEX = startsWith(string(rxnNames),'EX_');
rxn_isDM = startsWith(string(rxnNames),'DM_');
rxn_isTEX = contains(string(rxnNames),'tex');
rxn_isDrug = contains(string(rxnNames),'CM') | contains(string(rxnNames),'RFAMP') | contains(string(rxnNames),'DOXRBCN') | contains(string(rxnNames),'TTRCYC') | contains(string(rxnNames),'NOVBCN') | contains(string(rxnNames),'FUSA');
rxn_isBoundary = rxn_isEX | rxn_isDM | rxn_isTEX;
rxn_hasGPR = false(n,1);
if isfield(model_ir,'grRules') && numel(model_ir.grRules)==n
    rr = model_ir.grRules(:);
    rxn_hasGPR = cellfun(@(x) ~isempty(x) && any(isstrprop(x,'alphanum')), rr);
elseif isfield(model_ir,'rxnGeneMat') && isfield(model_ir,'genes')
    rxn_hasGPR = any(model_ir.rxnGeneMat~=0,2);
end
rxn_isInternal = ~rxn_isBoundary;

T_summary = table((1:n).', string(rxnNames(:)), rxn_hits, rxn_freq, rxn_isInternal, rxn_isBoundary, rxn_isEX, rxn_isDM, rxn_isTEX, rxn_isDrug, rxn_hasGPR, ...
    'VariableNames',{'rxnIdx','rxn','hits','freq','isInternal','isBoundary','isEX','isDM','isTEX','isDrug','hasGPR'});
writetable(T_summary, fullfile(outdir,'reaction_coverage_summary.csv'));

covered_any = rxn_hits>0;
miss_idx = find(~covered_any);
writematrix(string(rxnNames(miss_idx)), fullfile(outdir,'reactions_uncovered.txt'));

if exist('EFM_support','var') && ~isempty(EFM_support)
    supp_bin = logical(EFM_support);
else
    supp_bin = active;
end

badPartner = zeros(n,1);
for i = 1:n
    rn = rxnNames{i};
    if endsWith(rn,'_f')
        j = find(strcmp(rxnNames,[rn(1:end-2) '_b']),1);
        if ~isempty(j)
            badPartner(i)=j; badPartner(j)=i;
        end
    end
end
bpairs = find(badPartner>0 & (1:n)' < badPartner);
p = numel(bpairs);

pair_active = false(p,numEFMs);
pair_names = strings(p,1);
for t=1:p
    i = bpairs(t); j = badPartner(i);
    pair_active(t,:) = supp_bin(i,:) | supp_bin(j,:);
    pair_names(t) = sprintf('%s OR %s', rxnNames{i}, rxnNames{j});
end

pair_hits = sum(pair_active,2);
pair_freq = pair_hits/numEFMs;
irrevs = find(badPartner==0);
numIr = numel(irrevs);

group_names = [pair_names; string(rxnNames(irrevs))];
group_active = [pair_active; supp_bin(irrevs,:)];
group_hits = sum(group_active,2);
group_freq = group_hits/numEFMs;
group_covered = group_hits>0;

T_pairs = table(pair_names, pair_hits, pair_freq, 'VariableNames',{'pair','hits','freq'});
writetable(T_pairs, fullfile(outdir,'reversible_pair_coverage.csv'));

T_groups = table(group_names, group_hits, group_freq, group_covered, 'VariableNames',{'group','hits','freq','covered'});
writetable(T_groups, fullfile(outdir,'collapsed_group_coverage.csv'));

f = figure('Name','EFM length histogram'); histogram(efm_len); xlabel('Active reactions per EFM'); ylabel('Count'); title(sprintf('EFM lengths (n=%d)',numEFMs));
saveas(f, fullfile(outdir,'hist_efm_length.png')); close(f);

f = figure('Name','Residual histogram'); histogram(residual); xlabel('||S v||_{\infty}'); ylabel('Count'); title('Steady-state residuals across EFMs');
saveas(f, fullfile(outdir,'hist_residuals.png')); close(f);

f = figure('Name','Reaction coverage histogram'); histogram(rxn_hits); xlabel('EFM count per reaction'); ylabel('Number of reactions'); title('Reaction coverage');
saveas(f, fullfile(outdir,'hist_rxn_coverage.png')); close(f);

f = figure('Name','Collapsed group coverage histogram'); histogram(group_hits); xlabel('EFM count per group'); ylabel('Number of groups'); title('Collapsed (rev pairs) coverage');
saveas(f, fullfile(outdir,'hist_group_coverage.png')); close(f);

f = figure('Name','Length vs residual'); scatter(efm_len, residual, 'filled'); xlabel('Active reactions'); ylabel('||S v||_{\infty}'); title('EFM length vs residual');
saveas(f, fullfile(outdir,'scatter_length_vs_residual.png')); close(f);

[~,ix_top_rxn] = sort(rxn_hits,'descend');
kshow = min(30, n);
T_top_rxn = table(string(rxnNames(ix_top_rxn(1:kshow))), rxn_hits(ix_top_rxn(1:kshow)), rxn_freq(ix_top_rxn(1:kshow)), ...
    'VariableNames',{'rxn','hits','freq'});
writetable(T_top_rxn, fullfile(outdir,'top30_reactions_by_hits.csv'));

[~,ix_top_pair] = sort(pair_hits,'descend');
kshowp = min(30, p);
T_top_pair = table(pair_names(ix_top_pair(1:kshowp)), pair_hits(ix_top_pair(1:kshowp)), pair_freq(ix_top_pair(1:kshowp)), ...
    'VariableNames',{'pair','hits','freq'});
writetable(T_top_pair, fullfile(outdir,'top30_pairs_by_hits.csv'));

if isfield(model_ir,'subSystems') && numel(model_ir.subSystems)==n
    rawSubs = model_ir.subSystems(:);
    subs = strings(n,1);
    for k = 1:n
        x = rawSubs{k};
        if isempty(x)
            subs(k) = "NA";
        elseif iscell(x) || isstring(x)
            subs(k) = strjoin(string(x), '; ');
        elseif ischar(x)
            subs(k) = string(x);
        else
            subs(k) = string(x);
        end
    end
else
    subs = repmat("NA",n,1);
end

cats = ["Internal","Boundary","EX","DM","tex","Drug","GPR","NoGPR"];
mask.Internal = rxn_isInternal;
mask.Boundary = rxn_isBoundary;
mask.EX = rxn_isEX;
mask.DM = rxn_isDM;
mask.tex = rxn_isTEX;
mask.Drug = rxn_isDrug;
mask.GPR = rxn_hasGPR;
mask.NoGPR = ~rxn_hasGPR;

stats_cat = strings(0,1); stats_n = []; stats_cov = []; stats_frac = [];
fn = fieldnames(mask);
for u = 1:numel(fn)
    M = mask.(fn{u});
    stats_cat(end+1,1) = string(fn{u});
    stats_n(end+1,1) = sum(M);
    stats_cov(end+1,1) = sum(rxn_hits(M)>0);
    stats_frac(end+1,1) = stats_cov(end)/max(1,stats_n(end));
end
T_stats = table(stats_cat, stats_n, stats_cov, stats_frac, 'VariableNames',{'category','n','covered','covered_frac'});
writetable(T_stats, fullfile(outdir,'coverage_by_category.csv'));

if any(subs~="NA")
    us = unique(subs);
    subs_n = zeros(numel(us),1);
    subs_cov = zeros(numel(us),1);
    subs_hits_mean = zeros(numel(us),1);
    for sidx=1:numel(us)
        m = subs==us(sidx);
        subs_n(sidx) = sum(m);
        subs_cov(sidx) = sum(rxn_hits(m)>0);
        tmp = rxn_hits(m);
        if isempty(tmp), subs_hits_mean(sidx)=0; else subs_hits_mean(sidx)=mean(tmp); end
    end
    T_subs = table(us, subs_n, subs_cov, subs_cov./max(1,subs_n), subs_hits_mean, 'VariableNames',{'subsystem','n','covered','covered_frac','mean_hits'});
    writetable(T_subs, fullfile(outdir,'coverage_by_subsystem.csv'));
end

[~,sort_rxn] = sort(rxn_hits,'descend');
[~,sort_efm] = sort(efm_len,'descend');
Ashow = active(sort_rxn, sort_efm);
f = figure('Name','Heatmap EFMs x Reactions'); imagesc(Ashow); axis tight; xlabel('EFMs (sorted by length)'); ylabel('Reactions (sorted by frequency)'); title('Activity heatmap'); colorbar;
saveas(f, fullfile(outdir,'heatmap_active.png')); close(f);

cooc = active*active.';
diagv = diag(cooc);
cooc_norm = cooc;
for i=1:n
    if diagv(i)>0, cooc_norm(i,:) = cooc(i,:)/diagv(i); end
end
f = figure('Name','Co-occurrence heatmap (top 100 reactions)'); topN = min(100,n); idxTop = sort_rxn(1:topN); imagesc(cooc_norm(idxTop,idxTop)); axis square; colorbar; title('Reaction co-occurrence (row-normalized)'); xticks(1:topN); yticks(1:topN); xticklabels(rxnNames(idxTop)); yticklabels(rxnNames(idxTop)); xtickangle(90);
saveas(f, fullfile(outdir,'heatmap_cooccurrence_top100.png')); close(f);

maxJ = 500;
if numEFMs<=maxJ
    J = zeros(numEFMs);
    for i=1:numEFMs
        ai = active(:,i);
        for j=i:numEFMs
            aj = active(:,j);
            u = nnz(ai|aj); if u==0, J(i,j)=1; else, J(i,j)=nnz(ai&aj)/u; end
            J(j,i)=J(i,j);
        end
    end
    f = figure('Name','Jaccard EFMs'); imagesc(J); axis square; colorbar; title('Jaccard similarity between EFMs'); xlabel('EFM'); ylabel('EFM');
    saveas(f, fullfile(outdir,'heatmap_jaccard.png')); close(f);
    writematrix(J, fullfile(outdir,'jaccard_matrix.csv'));
end

if exist('efmAnchors','var') && ~isempty(efmAnchors)
    anchors = efmAnchors(:);
    anames = string(rxnNames(anchors));
    [ua,~,ic] = unique(anames);
    acount = accumarray(ic,1);
    T_anchor = table(ua, acount, 'VariableNames',{'anchor_rxn','count'});
    T_anchor = sortrows(T_anchor,'count','descend');
    writetable(T_anchor, fullfile(outdir,'anchor_reaction_counts.csv'));
end

badpairs_covered = any(pair_active,2);
bp_miss = find(~badpairs_covered);
if ~isempty(bp_miss)
    fid = fopen(fullfile(outdir,'reversible_pairs_uncovered.txt'),'w');
    for t=reshape(bp_miss,1,[])
        fprintf(fid,'%s\n', pair_names(t));
    end
    fclose(fid);
end

E_stats = table((1:numEFMs).', efm_len, residual, 'VariableNames',{'EFM','length','residual'});
writetable(E_stats, fullfile(outdir,'efm_stats.csv'));

save(fullfile(outdir,'analysis_workspace.mat'), 'rxn_hits','rxn_freq','efm_len','residual','active','pair_active','group_active','badPartner','bpairs','pair_names','group_names','cooc','cooc_norm');

disp('Analysis complete. Outputs in folder:'); disp(outdir);
end
