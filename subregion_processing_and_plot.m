%% Subregion analyses
%used for analyses in Figure 4 and Supplementary Figures related to
%partitioning of caudoputamen into four cortical input-defined clusters
%based on analysis of Hunnicutt et al. (2016). Relies on code from Tianyi
%Mao lab available on their GitHub under:
%https://github.com/TianyiMaoLab/striatumclusters. 
%Input data files available on Zenodo: https://zenodo.org/records/18685812

%% 0. Setup and Parameters
addpath('C:\Users\walki\Box\grp-psom-fuccillo-lab\Interneuron Anatomy\Scripts');
workingDir  = 'C:\Users\walki\Box\grp-psom-fuccillo-lab\Interneuron Anatomy\Brain Analysis\Brains';
cd(workingDir);

REGION = 'Caudoputamen';
ctKeys = {'TH','SST','PV'};                 % cell types to include

% Coordinate columns
COL_ML  = 8;
COL_DV  = 9;
COL_AP  = 10;

% Midline in CCFv3 (µm)
midline = 5700;

% Filename patterns:
% raw:            TH_04_F_points           or  TH_04_F_new_points
pat_raw = '^(?<ct>[A-Za-z\+]+)_(?<id>[0-9A-Za-z]+)_(?<sex>[MFmf])(?:_new)?_points$';

% classified:     PV_10222_F_points_classified_4clusters
%             or  TH_F02_new_points_classified_4clusters
pat_classified = [ ...
    '^(?<ct>[A-Za-z\+]+)_' ...
    '(?<id>[0-9A-Za-z]+)_' ...
    '(?<sex>[MFmf])' ...
    '(?:_new)?_' ...
    'points_classified_(?<k>[0-9]+)clusters$' ...
    ];

%% 1. Set up four subregions (dorsomedial caudoputamen, lateral caudoputamen, tail caudoputamen, anterior ventrolateral caudoputamen)
% Setup for subregions

subClusters  = [5 7 12 15];                 % IDs in the 4-cluster volume
clusterNames = {'DMS','DLS','Tail','VS'};   % labels
nSub = numel(subClusters);

% 4-cluster volume path (25 µm, IDs 5/7/12/15)
clusterVolFile = 'C:\Users\walki\Box\grp-psom-fuccillo-lab\Interneuron Anatomy\Scripts\striatum_04_clusters_ccf.tif';
cpMaskFile = 'structure_672.nrrd';  % 0 = outside, 1 = inside

assert(~isempty(cpMaskFile), 'structure_672.nrrd not found on path');


% Load 4-cluster label volume (CCFv3, 25 µm)
info = imfinfo(clusterVolFile);
nZ   = numel(info);
nX   = info(1).Width;
nY   = info(1).Height;

clusterVol = zeros(nY, nX, nZ, 'uint8');   % (rows=Y, cols=X), uint8 just 8 bit integer
for z = 1:nZ
    clusterVol(:,:,z) = imread(clusterVolFile, z);
end

V = nrrdread(cpMaskFile);          % V: [z y x] = [AP DV ML]
cpMask = permute(V > 0, [2 3 1]);  % -> [y x z] = [DV ML AP], matches clusterVol


% Zero out voxels outside CP
clusterVol(~cpMask) = 0;

voxelSize25_um = 25;
voxelVol25_mm3 = (voxelSize25_um/1000)^3;

% Volume (mm^3) of each target cluster ID
clusterVol_mm3 = zeros(1, nSub);
for i = 1:nSub
    k    = subClusters(i);
    % k = 5, then 7, then 12, then 15
    %clusterVol(y,x,z) = subregion ID for a given voxel
    nvox = nnz(clusterVol == k); %nnz = number of nonzero, counts all the ones, we get a one when the clusterVol == subregion (True, 1)
    clusterVol_mm3(i) = nvox * voxelVol25_mm3;
    fprintf('%s (ID %d): %d voxels, %.4f mm^3\n', ...
        clusterNames{i}, k, nvox, clusterVol_mm3(i));
end

%% 2. Import cell subcluster assignment data
% Per-hemisphere, per-cell-type counts in each cluster (inside_mask only)

csvFiles = dir(fullfile(workingDir, '*_4clusters.xlsx'));

% hemiDataByType.(ct) = struct array with one entry per hemisphere sample
% fields: file, brainID, sex, hemi, counts (1 x nSub), density (1 x nSub)
hemiDataByType = struct();
for c = 1:numel(ctKeys)
    hemiDataByType.(ctKeys{c}) = struct( ...
        'file',    {}, ...
        'brainID', {}, ...
        'sex',     {}, ...
        'hemi',    {}, ...
        'counts',  {}, ...
        'density', {} );
end

for f = csvFiles'
    [~, base, ~] = fileparts(f.name);

    % Parse filename (raw or classified pattern)
    m = regexp(base, pat_raw, 'names', 'once');
    if isempty(m)
        m = regexp(base, pat_classified, 'names', 'once');
    end
    if isempty(m)
        continue
    end

    ct = normCellType(m.ct);
    if isempty(ct) || ~ismember(ct, ctKeys)
        continue
    end

    brainID = m.id;
    sex     = upper(m.sex);

    T = readtable(fullfile(f.folder, f.name));

    % Need Region, cluster_id, inside_mask
    if ~ismember('Region',      T.Properties.VariableNames) || ...
            ~ismember('cluster_id',  T.Properties.VariableNames) || ...
            ~ismember('inside_mask', T.Properties.VariableNames)
        continue
    end
    %
    % Filter: Caudoputamen & inside_mask == true
    %insideMask = T.inside_mask;
    %if ~islogical(insideMask)
    %    insideMask = logical(insideMask);
    %end

    keep = strcmpi(T.Region, REGION); % & insideMask;
    T    = T(keep, :);
    if isempty(T)
        continue
    end

    % cluster_id to numeric if it is not already
    cid = T.cluster_id;
    if ~isnumeric(cid)
        cid = str2double(string(cid));
    end
    T.cluster_id = cid;

    % --- Hemisphere split via ML coordinate (unflipped) ---
    TL = T(T{:,COL_ML} <  midline, :);   % left hemi
    TR = T(T{:,COL_ML} >= midline, :);   % right hemi

    % Left hemisphere entry
    if ~isempty(TL)
        countsL = zeros(1, nSub);
        for j = 1:nSub
            countsL(j) = sum(TL.cluster_id == subClusters(j));
        end
        hemiDataByType.(ct)(end+1) = struct( ...
            'file',    fullfile(f.folder, f.name), ...
            'brainID', brainID, ...
            'sex',     sex, ...
            'hemi',    'L', ...
            'counts',  countsL, ...
            'density', countsL ./ clusterVol_mm3 ); % One hemisphere volume
    end

    % Right hemisphere entry
    if ~isempty(TR)
        countsR = zeros(1, nSub);
        for j = 1:nSub
            countsR(j) = sum(TR.cluster_id == subClusters(j));
        end
        hemiDataByType.(ct)(end+1) = struct( ...
            'file',    fullfile(f.folder, f.name), ...
            'brainID', brainID, ...
            'sex',     sex, ...
            'hemi',    'R', ...
            'counts',  countsR, ...
            'density', countsR ./ clusterVol_mm3 );
    end
end


%% 3. Boxplots: per-hemisphere densities (cells/mm^3) by CELL TYPE (each tile), grouped by subregion (DMCP, LCP, TCP, AVMCP)

% Flag: set true to split by sex (F/M) within each cell type
splitBySex = true;

ylims = [0 3200];
clusterNames = ["DMS","LS","TS","aVMS"];
ctKeys = {'SST','PV','TH'};
if splitBySex, tiledlayout(1,6); end
for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    hemiStruct = hemiDataByType.(ct);

    if isempty(hemiStruct)
        fprintf('%s: no contributing hemispheres; skipping boxplot.\n', ct);
        continue
    end

    % Stack [nHemi x nSub] density matrix
    D = vertcat(hemiStruct.density);

    if ~splitBySex
        % --- pooled sexes ---
        figure('Name', sprintf('%s density by subregion (pooled sexes)', ct));
        boxplot(D, 'Labels', clusterNames);
        ylim(ylims);
        ylabel('Cells / mm^3');
        title(sprintf('%s: per-hemisphere densities in DMS(5), DLS(7), Tail(12), VS(15)', ct));

    else
        % --- split by sex ---
        sexVec = string({hemiStruct.sex}).';   % nHemi x 1, 'F'/'M'



        sexes = ["F","M"];
        for si = 1:2
            sx = sexes(si);
            keep = sexVec == sx;

            nexttile; hold on; box on;

            if ~any(keep)
                title(sprintf('%s — Sex %s (n=0)', ct, sx));
                ylim(ylims);
                set(gca, 'XTick', 1:numel(clusterNames), 'XTickLabel', clusterNames);
                ylabel('Cells / mm^3');
                continue
            end

            boxplot(D(keep,:), 'Labels', clusterNames);
            ylim(ylims);
            ylabel('Cells / mm^3');
            title(sprintf('%s — Sex %s (n=%d)', ct, sx, sum(keep)));
        end
    end
end


%% 4. Boxplots (FLIPPED): per-hemisphere densities by SUBREGION (each tile), grouped by CT (SST/PV/TH)
% Tiles: 2 x 4
%   Row 1: Females
%   Row 2: Males
%   Cols: DMS, LS, TS, aVMS
%
% Within each tile: CT groups (SST/PV/TH) as boxplots

ylims = [0 3200];
clusterNames = ["DMS","LS","TS","aVMS"];
ctKeys = {"SST","PV","TH"};

sexes = ["F","M"];                 % row order (top F, bottom M)
sexRowName = containers.Map(["F","M"], ["Females","Males"]);  %#ok<NASGU>

figure;
t = tiledlayout(2,4, 'TileSpacing','compact', 'Padding','compact');

for si = 1:numel(sexes)
    sx = sexes(si);

    for r = 1:numel(clusterNames)
        regionName = clusterNames(r);

        ax = nexttile; hold(ax,'on'); box(ax,'on');

        % Build pooled vectors: one column per CT, padded with NaNs to same length
        dataCols = cell(1, numel(ctKeys));
        nPerCT   = zeros(1, numel(ctKeys));

        for c = 1:numel(ctKeys)
            ct = ctKeys{c};
            H = hemiDataByType.(ct);

            if isempty(H)
                dataCols{c} = NaN(0,1);
                nPerCT(c) = 0;
                continue;
            end

            sexVec = string({H.sex}).';
            keep = (sexVec == sx);

            if ~any(keep)
                dataCols{c} = NaN(0,1);
                nPerCT(c) = 0;
                continue;
            end

            D = vertcat(H.density);         % nHemi x 4
            v = D(keep, r);                 % nKeep x 1 for this region
            v = v(:);

            dataCols{c} = v;
            nPerCT(c) = numel(v);
        end

        % Pad to same length for boxplot matrix
        maxN = max(cellfun(@numel, dataCols));
        if maxN == 0
            % empty tile: still format axes
            ylim(ax, ylims);
            set(ax, 'XTick', 1:numel(ctKeys), 'XTickLabel', string(ctKeys));
            ylabel(ax, 'Cells / mm^3');
        else
            M = NaN(maxN, numel(ctKeys));
            for c = 1:numel(ctKeys)
                v = dataCols{c};
                if ~isempty(v)
                    M(1:numel(v), c) = v;
                end
            end

            boxplot(ax, M, 'Labels', string(ctKeys));
            ylim(ax, ylims);
            ylabel(ax, 'Cells / mm^3');
        end

        % Titles: "DMS (Females)" etc.
        if sx == "F"
            sexLabel = "Females";
        else
            sexLabel = "Males";
        end
        title(ax, sprintf('%s (%s)', regionName, sexLabel));

        % Optional: print n per CT in top-left
        % nStr = sprintf('n: SST=%d  PV=%d  TH=%d', nPerCT(1), nPerCT(2), nPerCT(3));
        % text(ax, 0.02, 0.95, nStr, 'Units','normalized', ...
        %     'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        %     'FontSize',9);

    end
end


%% 5. ANOVA across subregion densities (per CT) with LME + Tukey post-hocs
allOmnibusP = nan(numel(ctKeys),1);              % collect omnibus p (Subregion) per CT
allCTnames  = strings(numel(ctKeys),1);
p_sub_all = nan(numel(ctKeys),1);
p_int_all = nan(numel(ctKeys),1);
ct_all    = strings(numel(ctKeys),1);


for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    S  = hemiDataByType.(ct);
    if isempty(S), fprintf('%s: no hemispheres; skipping ANOVA.\n', ct); continue; end
    allCTnames(c) = ct;

    % ===== Long table build (robust to density being 1×4 OR N×4 per hemisphere) =====

    rows = cell(0,4);  % {Density, Subregion, Sex, HemiID}

    for h = 1:numel(S)
        dens = S(h).density;     % expected either 1×4 or N×4
        if isempty(dens), continue; end

        % If dens came in as 4×N, transpose it
        if size(dens,2) ~= numel(clusterNames) && size(dens,1) == numel(clusterNames)
            dens = dens.';
        end

        % Now require N×4
        if size(dens,2) ~= numel(clusterNames)
            error('Unexpected density shape for %s_%s: got %dx%d, expected N×%d', ...
                S(h).brainID, S(h).hemi, size(dens,1), size(dens,2), numel(clusterNames));
        end

        nRep   = size(dens,1);
        hemiID = sprintf('%s_%s', S(h).brainID, S(h).hemi);
        sex    = char(S(h).sex);   % 'M' or 'F'

        % Append rows: one per (replicate row × subregion)
        for r = 1:nRep
            for sidx = 1:numel(clusterNames)
                rows(end+1,:) = {dens(r,sidx), clusterNames{sidx}, sex, hemiID}; %#ok<AGROW>
            end
        end
    end

    T = cell2table(rows, 'VariableNames', {'Density','Subregion','Sex','HemiID'});

    % Drop NaNs
    T = T(isfinite(T.Density), :);

    % Stable categories
    T.Subregion = categorical(T.Subregion, clusterNames);
    T.Sex       = categorical(T.Sex, {'F','M'});
    T.HemiID    = categorical(T.HemiID);

    % Interaction only if estimable (no empty Subregion×Sex cells)
    ctab = varfun(@numel, T, ...
    'InputVariables','Density', ...
    'GroupingVariables',{'Subregion','Sex'});

    hasEmptyCell = any(ctab.numel_Density == 0);


    if hasEmptyCell
        mdl = fitlme(T, 'Density ~ 1 + Subregion + Sex + (1|HemiID)');
        fprintf('[%s] Missing Subregion×Sex cells -> additive model used.\n', ct);
    else
        mdl = fitlme(T, 'Density ~ 1 + Subregion*Sex + (1|HemiID)');
    end


    A = anova(mdl);  % omnibus for Subregion

    % Robust omnibus extraction
    [p_sub, stat_sub]   = get_anova_term_p(A, "Subregion");
    [p_sex, stat_sex]   = get_anova_term_p(A, "Sex");
    [p_int, stat_int]   = get_anova_term_p(A, "Subregion:Sex");
    fprintf('['); fprintf(ct); fprintf(']\n');
    fprintf('Subregion (LME):      stat=%s, p=%s\n', stat_sub, ternary(isfinite(p_sub), sprintf('%.3g',p_sub), 'NA'));
    fprintf('Sex (LME):           stat=%s, p=%s\n', stat_sex, ternary(isfinite(p_sex), sprintf('%.3g',p_sex), 'NA'));
    fprintf('Subregion×Sex (LME): stat=%s, p=%s\n', stat_int, ternary(isfinite(p_int), sprintf('%.3g',p_int), 'NA'));

    ct_all(c)    = ct;
    p_sub_all(c) = p_sub;
    p_int_all(c) = p_int;

    % Normalize accessors for table/dataset
    if istable(A)
        termNames = string(A.Properties.RowNames);
        varNames  = string(A.Properties.VariableNames);
        getcol = @(name) A.(name);
    else
        % dataset/titled dataset
        varNames  = string(get(A,'VarNames'));
        firstCol  = varNames(1);                    % first column holds term labels
        termNames = string(A.(firstCol));           % use dot, NOT braces
        getcol = @(name) A.(char(name));
    end

    % Find 'Subregion' row
    rowSub = strcmpi(strtrim(termNames), "Subregion");
    if ~any(rowSub)
        warning('ANOVA output has no ''Subregion'' row; setting omnibus p = NaN.');
        F_sub = NaN; p_sub = NaN;
    else
        % Locate p column(s) across versions
        pCandidates = lower(["pValue","p","ProbF","Prob_F","pvalue","Pr>F","Pr_Greater_F","ProbChiSq","Prob_Chi_Sq"]);
        pVar = "";
        for nm = varNames
            if any(strcmpi(lower(nm), pCandidates))
                pVar = nm; break
            end
        end
        if pVar == ""
            warning('Could not find p-value column in ANOVA table; setting p = NaN.');
            p_sub = NaN;
        else
            tmp = getcol(pVar);
            p_sub = double(tmp(rowSub));
            if numel(p_sub) ~= 1, p_sub = p_sub(1); end
        end

        % Try to get F or ChiSq for logging (optional)
        fCandidates = lower(["F","FStat","FStatistic","ChiSq","Chi_Square","ChiSquare"]);
        fVar = "";
        for nm = varNames
            if any(strcmpi(lower(nm), fCandidates))
                fVar = nm; break
            end
        end
        if fVar ~= ""
            tmp = getcol(fVar);
            F_sub = double(tmp(rowSub));
            if numel(F_sub) ~= 1, F_sub = F_sub(1); end
        else
            F_sub = NaN;
        end
    end

    fprintf('Subregion (LME): F/ChiSq = %s, p = %s\n', ...
        ternary(isfinite(F_sub), sprintf('%.3g',F_sub), 'NA'), ...
        ternary(isfinite(p_sub), sprintf('%.3g',p_sub), 'NA'));

    % Safe assignment (won’t error if p_sub is empty)
    if ~isempty(p_sub) && isscalar(p_sub)
        allOmnibusP(c) = p_sub;
    else
        allOmnibusP(c) = NaN;
    end



    % keep for cross-CT FDR later
    allOmnibusP(c) = p_sub;

    % =========================
    % Post-hoc pairwise (LME) via contrasts with coefTest
    % =========================
    subCats = categories(sub);
    K       = numel(subCats);

    % Fixed effects and their covariance
    beta   = fixedEffects(mdl);                % column vector
    FEcov  = mdl.CoefficientCovariance;        % covariance of fixed effects
    feNames = string(mdl.CoefficientNames);    % e.g., "(Intercept)","Subregion_A","Subregion_B",...

    % Identify which subregion is the reference (has no explicit coefficient)
    % Convention: treatment coding, first level is reference unless changed.
    % Build a map: level -> coefficient index (0 if reference)
    level2idx = containers.Map('KeyType','char','ValueType','int32');
    refLevel  = subCats{1};   % assumes current ordering; reordercats() earlier if needed
    for k = 1:K
        lvl = subCats{k};
        label = "Subregion_" + string(lvl);
        hit = find(feNames == label, 1);
        if isempty(hit), level2idx(char(lvl)) = 0; else, level2idx(char(lvl)) = int32(hit); end
    end

    % All pairwise contrasts (estimate, CI, Tukey-ish Sidak p-adjust)
    pairs  = nchoosek(1:K,2);
    npairs = size(pairs,1);
    alpha  = 0.05;
    alpha_sidak = 1 - (1 - alpha)^(1/npairs);   % FWER control within this CT
    df     = mdl.DFE;

    fprintf('[%s] LME pairwise (Subregion) via contrasts (Sidak-adjusted):\n', ct);
    for t = 1:npairs
        i = pairs(t,1); j = pairs(t,2);
        idx_i = level2idx(char(subCats{i}));
        idx_j = level2idx(char(subCats{j}));

        % Contrast vector L over fixed effects (size = length(beta))
        L = zeros(size(beta));
        if idx_i > 0, L(idx_i) =  1; end
        if idx_j > 0, L(idx_j) = -1; end
        % (Intercept) drops out automatically.

        est = L'*beta;
        se  = sqrt(L'*FEcov*L);
        tstat = est / se;
        p_raw = 2*tcdf(-abs(tstat), df);
        p_sidak = 1 - (1 - p_raw)^npairs;

        tcrit = tinv(1 - alpha_sidak/2, df);
        ci_lo = est - tcrit*se;
        ci_hi = est + tcrit*se;

        fprintf('  %s vs %s: Δ=%.3g (95%% CI [%0.3g, %0.3g]), p_sidak=%.3g\n', ...
            subCats{i}, subCats{j}, est, ci_lo, ci_hi, p_sidak);
    end


    % --- Post-hoc pairwise comparisons (Tukey) on Subregion within the LME ---
    % Prefer the built-in multcompare for LinearMixedModel (R2022b+). Fallback to contrasts.
    didPosthoc = false;
    try
        % Newer MATLABs support: multcompare(mdl,'Subregion',...)
        PH = multcompare(mdl, 'Subregion', 'ComparisonType','tukey-kramer', 'Display','off');
        % PH is a table with Group1, Group2, Estimate, SE, DF, tStat, pValue, CI
        fprintf('[%s] LME Tukey–Kramer pairwise (Subregion):\n', ct);
        for i = 1:height(PH)
            fprintf('  %s vs %s: Δ=%.3g (95%% CI [%0.3g, %0.3g]), p = %.3g\n', ...
                string(PH.Group1(i)), string(PH.Group2(i)), ...
                PH.Estimate(i), PH.Lower(i), PH.Upper(i), PH.pValue(i));
        end
        didPosthoc = true;
    catch
        % Fallback: pairwise contrasts on fixed effects via coefTest (version-robust)
        warning('%s: multcompare on LME unavailable; using manual contrasts (coefTest).', ct);

        % Levels in Subregion *as used in the model*
        subCats = categories(sub);
        K = numel(subCats);

        % Fixed-effects names and coefficients
        beta  = fixedEffects(mdl);                 % column vector
        bname = string(mdl.CoefficientNames);      % includes '(Intercept)', 'Subregion_<level>', ...

        % Determine reference level for Subregion (the one without a coefficient term)
        subCoefMask  = startsWith(bname, "Subregion_");
        subCoefNames = erase(bname(subCoefMask), "Subregion_");
        refLevel     = setdiff(subCats, subCoefNames);  % the level not listed gets coded as reference
        if isempty(refLevel)
            % Fallback if setdiff fails (rare): assume first category is reference
            refLevel = subCats(1);
        else
            refLevel = refLevel(1);
        end

        % Precompute index lookups for each non-reference level
        idxByLevel = containers.Map;
        for k = 1:numel(subCats)
            lev = subCats{k};
            nm  = "Subregion_" + string(lev);
            idx = find(bname == nm, 1);
            if ~isempty(idx), idxByLevel(lev) = idx; end
        end

        % Familywise adjustment (Sidak) across all pairs
        npairs = nchoosek(K, 2);
        alpha  = 0.05;
        alpha_sidak = 1 - (1 - alpha)^(1/npairs);

        fprintf('[%s] LME pairwise (manual via coefTest, Sidak-corrected):\n', ct);

        t = 0;
        for i = 1:K-1
            for j = i+1:K
                t = t + 1;

                % Build contrast L so that L*beta = mu_i - mu_j
                % With reference coding: mu_level = b0 + b(Subregion_level) if level ~= ref, else b0
                L = zeros(1, numel(bname));
                li = subCats{i};
                lj = subCats{j};

                if isKey(idxByLevel, li), L(idxByLevel(li)) =  1; end   % +b_i
                if isKey(idxByLevel, lj), L(idxByLevel(lj)) = -1; end   % -b_j

                % Test contrast
                [p_raw, Fstat, df1, df2] = coefTest(mdl, L, 0);  %#ok<ASGLU> df1 should be 1
                est = L * beta;                   % difference of marginal means
                se  = abs(est) / sqrt(max(Fstat, eps));  % from F = (est^2)/(se^2)

                % Sidak-adjusted p and CI
                p_adj = 1 - (1 - p_raw)^npairs;
                tcrit = tinv(1 - alpha_sidak/2, df2);
                ci_lo = est - tcrit * se;
                ci_hi = est + tcrit * se;

                fprintf('  %s vs %s: Δ=%.3g (95%% CI [%0.3g, %0.3g]), p_sidak = %.3g (raw p = %.3g)\n', ...
                    li, lj, est, ci_lo, ci_hi, p_adj, p_raw);
            end
        end
        didPosthoc = true;
    end


    % (Optional) keep any per-CT posthoc tables if you want to aggregate later
end

% --- Cross-CT multiple-testing correction on LME omnibus p-values (recommended)
ok = isfinite(allOmnibusP);
if any(ok)
    pvals = allOmnibusP(ok);
    [adjP, qval] = fdr_bh_local(pvals); % Benjamini–Hochberg
    validNames = allCTnames(ok);

    fprintf('\n== Omnibus Subregion effect across CTs (BH-FDR) ==\n');
    for k = 1:numel(pvals)
        fprintf('  %s: raw p = %.3g, FDR q = %.3g\n', validNames{k}, pvals(k), qval(k));
    end
end

alphaFDR = 0.05;
idx = find(ok);
q_int = nan(size(p_int_all));
q_int(idx) = bh_fdr_qvalues(p_int_all(idx));
sig_int = q_int <= alphaFDR;

fprintf('\n== Subregion×Sex omnibus across CTs (BH-FDR q=%.2f) ==\n', alphaFDR);
for ii = 1:numel(idx)
    c = idx(ii);
    fprintf('  %s: raw p=%.3g, FDR q=%.3g, sig=%d\n', ct_all(c), p_int_all(c), q_int(c), sig_int(c));
end
idx = find(isfinite(p_sub_all));
q_sub = nan(size(p_sub_all));
q_sub(idx) = bh_fdr_qvalues(p_sub_all(idx));
sig_sub = q_sub <= alphaFDR;

fprintf('\n== Subregion omnibus across CTs (BH-FDR q=%.2f) ==\n', alphaFDR);
for ii = 1:numel(idx)
    c = idx(ii);
    fprintf('  %s: raw p=%.3g, FDR q=%.3g, sig=%d\n', ct_all(c), p_sub_all(c), q_sub(c), sig_sub(c));
end


%% 6. ANOVA across cell types (per Subregion) with LME + Tukey post-hocs
allRegionP = nan(numel(clusterNames),1);
allRegionNames = clusterNames;
p_ct_all = nan(numel(clusterNames),1);
p_sex_all = nan(numel(clusterNames),1);
p_int_all = nan(numel(clusterNames),1);
sub_all = strings(numel(clusterNames),1);

for s = 1:numel(clusterNames)
    region = clusterNames{s};
    fprintf('\n==== Subregion %s ====\n', region);

    % Gather densities for this subregion across all CTs
    y = []; ct = []; hemiID = []; sex = [];
    for c = 1:numel(ctKeys)
        ctname = ctKeys{c};
        S = hemiDataByType.(ctname);
        if isempty(S), continue; end

        densVec = vertcat(S.density);           % [nHemi × nSub]
        nHemi = size(densVec,1);

        % Select the column corresponding to this subregion
        y_c = densVec(:, s);                    % density for subregion s
        hemi_c = strings(nHemi,1);
        for h = 1:nHemi
            hemi_c(h) = sprintf('%s_%s', S(h).brainID, S(h).hemi);
        end
        sex_c = strings(nHemi,1);
        for h = 1:nHemi
            sex_c(h) = S(h).sex;   % 'M' or 'F'
        end


        y = [y; y_c];
        hemiID = [hemiID; hemi_c];
        sex = [sex; sex_c];
        ct = [ct; repmat({ctname}, nHemi, 1)];
    end

    % Clean up
    ok = isfinite(y);
    y = y(ok); ct = categorical(ct(ok)); hemiID = categorical(hemiID(ok));
    sex = categorical(sex(ok), {'F','M'});

    % Mixed-effects model
    T = table(y, ct, sex, hemiID, ...
    'VariableNames', {'Density','CellType','Sex','HemiID'});
    mdl = fitlme(T, ...
    'Density ~ 1 + CellType * Sex + (1|HemiID)', ...
    'DummyVarCoding','reference');
    A = anova(mdl);
    disp(A)
    [p_ct,  ~] = get_anova_term_p(A, "CellType");
    [p_sex, ~] = get_anova_term_p(A, "Sex");
    [p_int, ~] = get_anova_term_p(A, "CellType:Sex");

    sub_all(s)    = region;
    p_ct_all(s) = p_ct;
    p_sex_all(s) = p_sex;
    p_int_all(s) = p_int;


    % Extract omnibus p for CellType
    if isa(A,'dataset')
        vn = get(A,'VarNames');
        termVar = vn{1};
        terms = string(A.(termVar));
        rowCellType = strcmpi(strtrim(terms),"CellType");
        pVar = vn{find(ismember(lower(vn), {'pvalue','p','probf','prob_f'}),1)};
        p_sub = A.(pVar)(rowCellType);
    else
        rowCellType = strcmpi(A.Properties.RowNames,"CellType");
        p_sub = A.pValue(rowCellType);
    end
    allRegionP(s) = p_sub;
    fprintf('Subregion %s: omnibus CellType effect p = %.3g\n', region, p_sub);

    % --- Post-hoc pairwise between cell types (per subregion) ---
    ctCats = categories(ct);
    K = numel(ctCats);

    % Fixed effects and covariance (coerce to numeric if needed)
    beta = fixedEffects(mdl);                    % p×1
    CovB = mdl.CoefficientCovariance;
    if ~isnumeric(CovB)
        if istable(CovB),    CovB = table2array(CovB);
        elseif isa(CovB,'dataset') || isa(CovB,'titleddataset')
            CovB = table2array(dataset2table(CovB));
        else,                CovB = double(CovB);   % last resort
        end
    end

    coefNames = string(mdl.CoefficientNames);    % e.g., "(Intercept)","CellType_TH",...

    % Map each CellType level to its coefficient index (empty -> reference level)
    lvl2idx = containers.Map('KeyType','char','ValueType','any');
    for k = 1:K
        lvl = ctCats{k};
        % find name that looks like "CellType_<level>" (robust to spaces/encoding)
        hit = find(contains(coefNames, "CellType_", 'IgnoreCase',true) & ...
            contains(coefNames, lvl,        'IgnoreCase',true), 1);
        if ~isempty(hit)
            lvl2idx(char(lvl)) = hit;   % non-reference level has a coefficient
        else
            lvl2idx(char(lvl)) = [];    % reference level (no separate coefficient)
        end
    end

    alpha    = 0.05;
    npairs   = nchoosek(K,2);
    alphaSid = 1 - (1 - alpha)^(1/npairs);
    df       = mdl.DFE;

    fprintf('[%s] LME pairwise (CellType) via fixed-effect contrasts (Sidak):\n', region);
    for i = 1:K-1
        for j = i+1:K
            Li = zeros(1, numel(beta));
            Lj = zeros(1, numel(beta));

            idx_i = lvl2idx(char(ctCats{i}));
            idx_j = lvl2idx(char(ctCats{j}));
            if ~isempty(idx_i), Li(idx_i) =  1; end
            if ~isempty(idx_j), Lj(idx_j) =  1; end

            % Contrast for mean difference (i - j); intercept cancels automatically
            L   = Li - Lj;                      % 1×p
            est = L * beta;                     % scalar
            se  = sqrt(L * CovB * L.');         % scalar

            % If a level mapping failed for some reason, fall back to coefTest
            if ~(isfinite(se) && se>0)
                [p_raw, F, df1, df2] = coefTest(mdl, L); %#ok<ASGLU>
                p_sid = 1 - (1 - p_raw)^npairs;
                fprintf('  %s vs %s: Δ=%.3g, p_sidak = %.3g (raw p = %.3g)\n', ...
                    ctCats{i}, ctCats{j}, est, p_sid, p_raw);
                continue;
            end

            tstat  = est / se;
            p_raw  = 2 * tcdf(-abs(tstat), df);
            p_sid  = 1 - (1 - p_raw)^npairs;

            tcrit  = tinv(1 - alphaSid/2, df);
            ci_lo  = est - tcrit*se;
            ci_hi  = est + tcrit*se;

            fprintf('  %s vs %s: Δ=%.3g (95%% CI [%0.3g, %0.3g]), p_sidak = %.3g (raw p = %.3g)\n', ...
                ctCats{i}, ctCats{j}, est, ci_lo, ci_hi, p_sid, p_raw);
        end

    end
end

% Cross-region FDR correction
ok = isfinite(allRegionP);
if any(ok)
    pvals = allRegionP(ok);
    [~, qvals] = fdr_bh_local(pvals);
    fprintf('\n== CellType effect across subregions (BH-FDR) ==\n');
    for k = 1:numel(pvals)
        fprintf('  %s: raw p = %.3g, q = %.3g\n', allRegionNames{k}, pvals(k), qvals(k));
    end
end

alphaFDR = 0.05;
idx = find(ok);
q_int = nan(size(p_int_all));
q_int(idx) = bh_fdr_qvalues(p_int_all(idx));
sig_int = q_int <= alphaFDR;

fprintf('\n== Subregion×Sex omnibus across CTs (BH-FDR q=%.2f) ==\n', alphaFDR);
for ii = 1:numel(idx)
    c = idx(ii);
    fprintf('  %s: raw p=%.3g, FDR q=%.3g, sig=%d\n', sub_all(c), p_int_all(c), q_int(c), sig_int(c));
end
idx = find(isfinite(p_ct_all));
q_ct = nan(size(p_ct_all));
q_ct(idx) = bh_fdr_qvalues(p_ct_all(idx));
sig_ct = q_ct <= alphaFDR;

fprintf('\n== Subregion omnibus across CTs (BH-FDR q=%.2f) ==\n', alphaFDR);
for ii = 1:numel(idx)
    c = idx(ii);
    fprintf('  %s: raw p=%.3g, FDR q=%.3g, sig=%d\n', sub_all(c), p_ct_all(c), q_ct(c), sig_ct(c));
end



followupRegions = {'DLS','VS'};
alpha = 0.05;

% (optional) enforce a consistent CT order in prints + contrasts
% If you want: ctKeys = {'TH','SST','PV'};  % uncomment if desired

for r = 1:numel(followupRegions)
    region = followupRegions{r};
    fprintf('\n===== FOLLOW-UP: %s =====\n', region);

    regionIdx = find(strcmp(clusterNames, region), 1);
    assert(~isempty(regionIdx), 'Region %s not found in clusterNames.', region);

    for sx = {'F','M'}
        sexStr = sx{1};
        fprintf('--- Sex = %s ---\n', sexStr);

        y = []; ct = {}; hemiID = strings(0,1);

        % ---------- gather data ----------
        for c = 1:numel(ctKeys)
            ctname = ctKeys{c};
            S = hemiDataByType.(ctname);
            if isempty(S), continue; end

            densMat = vertcat(S.density);      % nHemi × nSub
            sexVec  = string({S.sex}).';       % nHemi × 1

            keep = sexVec == sexStr;
            if ~any(keep), continue; end

            y_c = densMat(keep, regionIdx);    % densities for this region+sex
            idx = find(keep);

            hemi_c = strings(numel(idx),1);
            for k = 1:numel(idx)
                h = idx(k);
                hemi_c(k) = sprintf('%s_%s', S(h).brainID, S(h).hemi);
            end

            y     = [y; y_c];
            ct    = [ct; repmat({ctname}, numel(y_c), 1)];
            hemiID = [hemiID; hemi_c];
        end

        % Drop NaNs / infs (important before both summaries + fitlme)
        ok = isfinite(y);
        y = y(ok);
        ct = ct(ok);
        hemiID = hemiID(ok);

        if isempty(y)
            fprintf('  (no data for %s %s)\n', region, sexStr);
            continue;
        end

        % ---------- RAW SUMMARIES (the clutch part) ----------
        % Prints: n, mean±SD, median[IQR] PER cell type, using the exact data in this follow-up
        ctCats = unique(ct, 'stable');  % preserves ctKeys order in printout
        fprintf('Raw densities (cells/mm^3), by CellType (region %s, sex %s):\n', region, sexStr);

        for k = 1:numel(ctCats)
            thisCT = ctCats{k};
            yy = y(strcmp(ct, thisCT));

            if isempty(yy)
                fprintf('  %s: n=0\n', thisCT);
                continue
            end

            fprintf('  %s: n=%d | mean=%.2f ± %.2f | median=%.2f [IQR %.2f–%.2f]\n', ...
                thisCT, numel(yy), mean(yy), std(yy), median(yy), prctile(yy,25), prctile(yy,75));
        end

        % Optional: print ordering by mean and by median (quick sanity)
        mu = nan(numel(ctCats),1);
        md = nan(numel(ctCats),1);
        for k = 1:numel(ctCats)
            yy = y(strcmp(ct, ctCats{k}));
            mu(k) = mean(yy);
            md(k) = median(yy);
        end
        [~, ordMu] = sort(mu, 'descend');
        [~, ordMd] = sort(md, 'descend');
        fprintf('  Order by mean:   %s\n', strjoin(ctCats(ordMu), ' > '));
        fprintf('  Order by median: %s\n', strjoin(ctCats(ordMd), ' > '));

        % ---------- MODEL ----------
        T = table(y, categorical(ct, ctCats), categorical(hemiID), ...
            'VariableNames', {'Density','CellType','HemiID'});

        mdl = fitlme(T, 'Density ~ CellType + (1|HemiID)', ...
                     'DummyVarCoding','reference');

        % omnibus
        A = anova(mdl);
        [p_ct, stat_ct] = get_anova_term_p(A, "CellType");
        fprintf('CellType omnibus: stat=%s, p=%.3g\n', stat_ct, p_ct);

        % post-hoc only if omnibus is meaningful
        if p_ct < alpha
            fprintf('Pairwise CellType contrasts (%s, %s):\n', region, sexStr);

            ctCats = categories(T.CellType);   % now categorical order
            beta   = fixedEffects(mdl);
            CovB   = mdl.CoefficientCovariance;
            df     = mdl.DFE;
            npairs = nchoosek(numel(ctCats),2);
            alphaSid = 1 - (1 - alpha)^(1/npairs);

            coefNames = string(mdl.CoefficientNames);
            lvl2idx = containers.Map('KeyType','char','ValueType','any');
            for k = 1:numel(ctCats)
                hit = find(contains(coefNames, "CellType_") & ...
                           contains(coefNames, string(ctCats{k})), 1);
                if isempty(hit), lvl2idx(char(ctCats{k})) = [];
                else,            lvl2idx(char(ctCats{k})) = hit;
                end
            end

            for i = 1:numel(ctCats)-1
                for j = i+1:numel(ctCats)
                    L = zeros(1, numel(beta));
                    if ~isempty(lvl2idx(char(ctCats{i}))), L(lvl2idx(char(ctCats{i}))) =  1; end
                    if ~isempty(lvl2idx(char(ctCats{j}))), L(lvl2idx(char(ctCats{j}))) = -1; end

                    est = L * beta;
                    se  = sqrt(L * CovB * L.');
                    t   = est / se;
                    p_raw = 2 * tcdf(-abs(t), df);
                    p_sid = 1 - (1 - p_raw)^npairs;

                    fprintf('  %s vs %s: Δ=%.2f, p_sidak=%.3g (raw=%.3g)\n', ...
                        string(ctCats{i}), string(ctCats{j}), est, p_sid, p_raw);
                end
            end
        end
    end
end

%% 7. Plot DLS and VS: per-hemisphere densities by CellType, split by Sex
regionsToPlot = {'DLS','VS'};

% set CT order explicitly (so plots are consistent)
ctOrder = {'TH','SST','PV'};

% set y-lims if you want consistency across plots
ylims = [0 3200];

for r = 1:numel(regionsToPlot)
    region = regionsToPlot{r};
    sIdx = find(strcmp(clusterNames, region), 1);
    assert(~isempty(sIdx), 'Region %s not found in clusterNames.', region);

    % Build long table: one row = one hemisphere observation for that region
    y = []; ct = {}; sex = {}; hemiID = {};

    for c = 1:numel(ctOrder)
        ctname = ctOrder{c};
        S = hemiDataByType.(ctname);
        if isempty(S), continue; end

        D = vertcat(S.density);              % nHemi x nSub
        y_c = D(:, sIdx);                    % nHemi x 1 densities for this region

        nH = numel(S);
        hemi_c = strings(nH,1);
        sex_c  = strings(nH,1);
        for h = 1:nH
            hemi_c(h) = sprintf('%s_%s', S(h).brainID, S(h).hemi);
            sex_c(h)  = string(S(h).sex);    % 'M' or 'F'
        end

        y    = [y; y_c];
        ct   = [ct; repmat({ctname}, nH, 1)];
        sex  = [sex; cellstr(sex_c)];
        hemiID = [hemiID; cellstr(hemi_c)];
    end

    T = table(y, categorical(ct, ctOrder), categorical(sex, {'F','M'}), categorical(hemiID), ...
        'VariableNames', {'Density','CellType','Sex','HemiID'});

    % drop NaNs
    T = T(isfinite(T.Density), :);

    % ---- Plot: 2 panels (F, M) ----
    figure('Name', sprintf('%s: density by cell type, split by sex', region));

    for si = 1:2
        thisSex = categorical({'F','M'}); %#ok<NASGU>
    end

    sexes = {'F','M'};
    for si = 1:2
        sx = sexes{si};
        Ts = T(T.Sex == sx, :);

        subplot(1,2,si); hold on; box on;
        boxplot(Ts.Density, Ts.CellType, 'Labels', ctOrder);

        ylim(ylims);
        ylabel('Cells / mm^3');
        title(sprintf('%s — Sex %s', region, sx));

        % optional: overlay points (nice for n~small)
        % (jittered scatter)
        g = double(Ts.CellType);
        xj = g + (rand(size(g)) - 0.5) * 0.15;
        plot(xj, Ts.Density, '.', 'MarkerSize', 10);
    end
end


%% 8. Quick visualization to check if clusterVol is 1 hemisphere or whole brain
% Looks like this is one hemisphere

midZ = round(nZ/2);
figure;
imagesc(cpMask(:,:,midZ)); axis image off; title('cpMask slice'); colormap(gray);

figure;
imagesc(clusterVol(:,:,midZ)); axis image off; title('clusterVol slice'); colormap(jet); colorbar;

% Overlay (clusterVol on top of cpMask)
figure;
imagesc(cpMask(:,:,midZ)); axis image off; colormap(gray); hold on;
h = imagesc(clusterVol(:,:,midZ));
set(h,'AlphaData', clusterVol(:,:,midZ) > 0);  % only show labeled voxels
title('Overlay: clusterVol (labels) on cpMask');
%figure('Name','4-cluster volume — mid coronal slice');
%imagesc(clusterVol(:,:,midZ));          % raw labels: 0,5,7,12,15,...
%axis image off;
%colormap(jet);
%colorbar;
%title(sprintf('clusterVol mid-slice (z = %d)', midZ));

% Print unique labels in this slice & full volume for sanity
%u_mid  = unique(clusterVol(:,:,midZ));
%u_full = unique(clusterVol(:));
%fprintf('Mid-slice labels: '); fprintf('%d ', u_mid); fprintf('\n');
%fprintf('All-volume labels: '); fprintf('%d ', u_full); fprintf('\n');

%% 9. QC: plot one AP plane of detections colored by 5/7/12/15

ctQC    = 'PV';        % choose CT to inspect (e.g. 'TH','SST','PV')
brainQC = '102';     % choose brain ID to inspect

qcFile = '';
for f = csvFiles'
    [~, base, ~] = fileparts(f.name);
    m = regexp(base, pat_raw, 'names', 'once');
    if isempty(m)
        m = regexp(base, pat_classified, 'names', 'once');
    end
    if isempty(m), continue; end

    if strcmpi(normCellType(m.ct), ctQC) && strcmp(m.id, brainQC)
        qcFile = fullfile(f.folder, f.name);
        break;
    end
end

if ~isempty(qcFile)
    Tqc = readtable(qcFile);

    % CPu + inside_mask only
    insideMask = Tqc.inside_mask;
    if ~islogical(insideMask), insideMask = logical(insideMask); end
    keep = strcmpi(Tqc.Region, REGION) & insideMask;
    Tqc  = Tqc(keep, :);
    if ~isempty(Tqc)
        % cluster_id numeric + keep 5/7/12/15
        cid = Tqc.cluster_id;
        if ~isnumeric(cid), cid = str2double(string(cid)); end
        use = ismember(cid, subClusters);
        Tqc = Tqc(use, :);
        cid = cid(use);

        if ~isempty(Tqc)
            % choose one AP slab (150 µm) with most points
            ap = double(Tqc{:, COL_AP});
            edges = min(ap):150:max(ap);
            if numel(edges) >= 2
                [N,~,bin] = histcounts(ap, edges);
                [~,imax]  = max(N);
                slab = ap >= edges(imax) & ap < edges(imax+1);

                ml       = double(Tqc{slab, COL_ML});
                dv       = double(Tqc{slab, COL_DV});
                cid_slab = cid(slab);

                figure('Name', sprintf('QC %s %s — one AP plane', ctQC, brainQC));
                hold on; axis equal ij;
                colors = lines(numel(subClusters));
                for j = 1:numel(subClusters)
                    k = subClusters(j);
                    idx = (cid_slab == k);
                    if any(idx)
                        scatter(ml(idx), dv(idx), 30, colors(j,:), 'filled', ...
                            'DisplayName', sprintf('%s (%d)', clusterNames{j}, k));
                    end
                end
                xlabel('ML (µm)'); ylabel('DV (µm)');
                title('Detections in one AP slab (inside\_mask), colored by 4-cluster ID');
                legend('Location','bestoutside');
            end
        end
    end
end


%% Helper Functions
function s = normCellType(raw)
r = upper(regexprep(raw,'\+',''));  % drop '+' if present
switch r
    case 'CHAT', s = 'ChAT';
    case {'TH','SST','PV'}, s = r;
    otherwise, s = '';  % excluded/unsupported
end
end

% --- Local BH helper (no toolbox dependency)
function [padj, q] = fdr_bh_local(p)
%FDR_BH_LOCAL Benjamini–Hochberg FDR adjustment.
%   [padj, q] = fdr_bh_local(p)
%   p   : vector or array of raw p-values (NaNs allowed)
%   padj: BH-adjusted p-values (a.k.a. q-values)
%   q   : identical to padj for convenience

sz  = size(p);
p   = p(:);                     % vectorize
m   = numel(p);

% Keep track of finite entries
ok  = isfinite(p);
padj = nan(m,1);

if any(ok)
    % Sort the finite p-values
    [ps, idx] = sort(p(ok));
    r = (1:sum(ok))';           % ranks 1..m_ok
    m_ok = numel(ps);

    % BH raw q's
    qraw = ps .* (m_ok ./ r);

    % Enforce monotonicity from the end: flip -> cummin -> flip back
    qmono = flipud( cummin( flipud(qraw) ) );

    % Cap at 1
    qmono = min(qmono, 1);

    % Scatter back to original positions
    tmp = nan(sum(ok),1);
    tmp(idx) = qmono;
    padj(ok) = tmp;
end

% Reshape and mirror to q
padj = reshape(padj, sz);
q    = padj;
end

function y = min_accumulate_flip(x)
% ensure monotone non-increasing from end to start
y = x;
for i = numel(x)-1:-1:1
    y(i) = min(y(i), y(i+1));
end
end
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function [pTerm, statStr] = get_anova_term_p(A, termName)
% Extract p-value + a stat for a given term from anova(mdl) output.
% Compatible with table and titleddataset.

pTerm  = NaN;
statStr = "NA";

if istable(A)
    termNames = string(A.Properties.RowNames);
    varNames  = string(A.Properties.VariableNames);
    getcol = @(name) A.(name);
else
    % titleddataset
    varNames = string(get(A,'VarNames'));
    firstCol = varNames(1);
    termNames = string(A.(char(firstCol)));
    getcol = @(name) A.(char(name));
end

row = strcmpi(strtrim(termNames), string(termName));
if ~any(row), return; end

% p column candidates
pCandidates = ["pValue","p","ProbF","Prob_F","pvalue","Pr>F","ProbChiSq","Prob_Chi_Sq"];
pVar = "";
for nm = varNames
    if any(strcmpi(nm, pCandidates))
        pVar = nm; break
    end
end
if pVar ~= ""
    tmp = getcol(pVar);
    pTerm = double(tmp(row));
    if numel(pTerm) ~= 1, pTerm = pTerm(1); end
end

% stat column candidates (optional)
statCandidates = ["F","FStat","FStatistic","ChiSq","ChiSquare","Chi_Square"];
sVar = "";
for nm = varNames
    if any(strcmpi(nm, statCandidates))
        sVar = nm; break
    end
end
if sVar ~= ""
    tmp = getcol(sVar);
    statVal = double(tmp(row));
    if numel(statVal) ~= 1, statVal = statVal(1); end
    statStr = sprintf('%.3g', statVal);
end
end

function q = bh_fdr_qvalues(p)
%BH_FDR_QVALUES  Benjamini–Hochberg FDR q-values (adjusted p-values).
%   q = bh_fdr_qvalues(p) returns q-values the same size as p.
%
%   NaNs are preserved. p-values must be in [0,1].

q = nan(size(p));
ok = isfinite(p);
p = p(ok);

if isempty(p), return; end

p = double(p(:));
m = numel(p);

% sort
[ps, order] = sort(p, 'ascend');

% BH adjusted p-values
qsorted = ps .* (m ./ (1:m)');

% enforce monotonicity (step-down)
qsorted = flipud(cummin(flipud(qsorted)));

% cap at 1
qsorted(qsorted > 1) = 1;

% unsort back into original ok positions
q_ok = nan(m,1);
q_ok(order) = qsorted;

q(ok) = reshape(q_ok, size(q(ok)));
end