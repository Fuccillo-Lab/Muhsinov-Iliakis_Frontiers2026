%% Subregion analyses - additional
%used for analyses in Supplementary Figures related to
%partitioning of caudoputamen into four cortical input-defined clusters
%based on analysis of Hunnicutt et al. (2016). Relies on code from Tianyi
%Mao lab available on their GitHub under:
%https://github.com/TianyiMaoLab/striatumclusters. 
%Input data files available on Zenodo: https://zenodo.org/records/18685812

%% 0. Setup and Parameters
workingDir  = 'C:\Users\Jonib\Box\grp-psom-fuccillo-lab\Interneuron Anatomy\Brain Analysis\Brains';
cd(workingDir);

REGION  = 'Caudoputamen';
ctKeys  = {'TH','SST','PV'};        % cell types to include
sexKeys = {'F','M'};                % NEW: sexes to consider

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

subClusters  = [5 7 12 15];                 % IDs in the 4-cluster volume
clusterNames = {'DMS','DLS','Tail','VS'};   % labels
nSub = numel(subClusters);

% 4-cluster volume path (25 µm, IDs 5/7/12/15)
clusterVolFile = 'C:\Users\Jonib\Box\grp-psom-fuccillo-lab\Interneuron Anatomy\Scripts\striatum_04_clusters_ccf.tif';
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

    % Filter: Caudoputamen & inside_mask == true
    %insideMask = T.inside_mask;
    %if ~islogical(insideMask)
    %    insideMask = logical(insideMask);
    %end

    keep = strcmpi(T.Region, REGION);
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


%% 3. Boxplots: per-hemisphere densities (cells/mm^3) by CELL TYPE (each tile), grouped by subregion (DMCP, LCP, TCP, AVMCP), split by sex

for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    hemiStruct = hemiDataByType.(ct);

    if isempty(hemiStruct)
        fprintf('%s: no contributing hemispheres; skipping boxplot.\n', ct);
        continue
    end

    for s = 1:numel(sexKeys)            % NEW: loop over F/M
        sexKey = sexKeys{s};
        idxSex = strcmp({hemiStruct.sex}, sexKey);
        if ~any(idxSex)
            fprintf('%s (%s): no contributing hemispheres; skipping boxplot.\n', ct, sexKey);
            continue
        end

        % Stack [nHemi x nSub] density matrix for this sex
        D = vertcat(hemiStruct(idxSex).density);

        figure('Name', sprintf('%s density by subregion (per hemisphere, inside mask) — sex %s', ct, sexKey));
        boxplot(D, 'Labels', clusterNames);
        ylim([0 3200]);
        ylabel('Cells / mm^3');
        title(sprintf('%s (sex %s): per-hemisphere densities in DMS(5), DLS(7), Tail(12), VS(15)', ct, sexKey));
    end
end

%% 4. ANOVA across subregion densities (per CT) with LME + Tukey post-hocs, split by sex

allOmnibusP = nan(numel(ctKeys), numel(sexKeys));   % NEW: p per CT × sex
allCTnames  = string(ctKeys(:));

for c = 1:numel(ctKeys)
    ct    = ctKeys{c};
    S_all = hemiDataByType.(ct);
    if isempty(S_all)
        fprintf('%s: no hemispheres at all; skipping ANOVA.\n', ct);
        continue;
    end

    for sSex = 1:numel(sexKeys)         % NEW: loop over sex for stats
        sexKey = sexKeys{sSex};
        S = S_all(strcmp({S_all.sex}, sexKey));
        if isempty(S)
            fprintf('%s (sex %s): no hemispheres; skipping ANOVA.\n', ct, sexKey);
            continue;
        end

        % [nHemi x nSub] matrix of densities
        D     = vertcat(S.density);
        nHemi = size(D,1);
        nSub  = size(D,2);

        % --- Long-form vectors ---
        y   = D(:);  % densities
        sub = categorical(reshape(repmat(clusterNames, nHemi, 1), [], 1));  % subregion labels
        hemiID = strings(nHemi*nSub,1);
        for h = 1:nHemi
            hemiID((h-1)*nSub+1:h*nSub) = sprintf('%s_%s', S(h).brainID, S(h).hemi);
        end
        hemiID = categorical(hemiID);

        % Drop NaNs
        ok = isfinite(y);
        y = y(ok); sub = sub(ok); hemiID = hemiID(ok);

        fprintf('\n==== %s (sex %s) ====\n', ct, sexKey);

        % --- Quick one-way ANOVA (unchanged) ---
        [p, ~, stats] = anova1(y, sub, 'off'); %#ok<ASGLU>
        fprintf('[%s, sex %s] One-way ANOVA across subregions: p = %.3g\n', ct, sexKey, p);

        % --- Mixed-effects ANOVA (accounts for repeated measures) ---
        T   = table(y, sub, hemiID, 'VariableNames', {'Density','Subregion','HemiID'});
        mdl = fitlme(T, 'Density ~ 1 + Subregion + (1|HemiID)');

        A = anova(mdl);  % omnibus for Subregion

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

        fprintf('Subregion (LME, sex %s): F/ChiSq = %s, p = %s\n', ...
            sexKey, ...
            ternary(isfinite(F_sub), sprintf('%.3g',F_sub), 'NA'), ...
            ternary(isfinite(p_sub), sprintf('%.3g',p_sub), 'NA'));

        % Safe assignment
        if ~isempty(p_sub) && isscalar(p_sub)
            allOmnibusP(c, sSex) = p_sub;          % NEW: store by CT × sex
        else
            allOmnibusP(c, sSex) = NaN;
        end

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

        fprintf('[%s, sex %s] LME pairwise (Subregion) via contrasts (Sidak-adjusted):\n', ct, sexKey);
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
        try
            PH = multcompare(mdl, 'Subregion', 'ComparisonType','tukey-kramer', 'Display','off');
            fprintf('[%s, sex %s] LME Tukey–Kramer pairwise (Subregion):\n', ct, sexKey);
            for i = 1:height(PH)
                fprintf('  %s vs %s: Δ=%.3g (95%% CI [%0.3g, %0.3g]), p = %.3g\n', ...
                    string(PH.Group1(i)), string(PH.Group2(i)), ...
                    PH.Estimate(i), PH.Lower(i), PH.Upper(i), PH.pValue(i));
            end
        catch
            warning('%s (sex %s): multcompare on LME unavailable; using manual contrasts (coefTest).', ct, sexKey);

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

            fprintf('[%s, sex %s] LME pairwise (manual via coefTest, Sidak-corrected):\n', ct, sexKey);

            t = 0;
            for i = 1:K-1
                for j = i+1:K
                    t = t + 1;

                    % Build contrast L so that L*beta = mu_i - mu_j
                    L = zeros(1, numel(bname));
                    li = subCats{i};
                    lj = subCats{j};

                    if isKey(idxByLevel, li), L(idxByLevel(li)) =  1; end
                    if isKey(idxByLevel, lj), L(idxByLevel(lj)) = -1; end

                    [p_raw, Fstat, df1, df2] = coefTest(mdl, L, 0); %#ok<ASGLU>
                    est = L * beta;
                    se  = abs(est) / sqrt(max(Fstat, eps));

                    p_sidak = 1 - (1 - p_raw)^npairs;
                    tcrit   = tinv(1 - alpha_sidak/2, df2);
                    ci_lo   = est - tcrit * se;
                    ci_hi   = est + tcrit * se;

                    fprintf('  %s vs %s: Δ=%.3g (95%% CI [%0.3g, %0.3g]), p_sidak = %.3g (raw p = %.3g)\n', ...
                        li, lj, est, ci_lo, ci_hi, p_sidak, p_raw);
                end
            end
        end
    end
end

% --- Cross-CT multiple-testing correction on LME omnibus p-values (per sex)
for sSex = 1:numel(sexKeys)
    sexKey = sexKeys{sSex};
    pcol = allOmnibusP(:, sSex);
    ok = isfinite(pcol);
    if ~any(ok), continue; end

    pvals = pcol(ok);
    [~, qval] = fdr_bh_local(pvals); % Benjamini–Hochberg
    validNames = allCTnames(ok);

    fprintf('\n== Omnibus Subregion effect across CTs (BH-FDR) — sex %s ==\n', sexKey);
    for k = 1:numel(pvals)
        fprintf('  %s: raw p = %.3g, FDR q = %.3g\n', validNames(k), pvals(k), qval(k));
    end
end


%% 5. ANOVA across cell types (per Subregion) with LME + Tukey post-hocs, split by sex

allRegionP     = nan(numel(clusterNames), numel(sexKeys));   % NEW: p per subregion × sex
allRegionNames = clusterNames;

for rIdx = 1:numel(clusterNames)
    region = clusterNames{rIdx};

    for sSex = 1:numel(sexKeys)
        sexKey = sexKeys{sSex};
        fprintf('\n==== Subregion %s (sex %s) ====\n', region, sexKey);

        % Gather densities for this subregion across all CTs, this sex only
        y = []; ct = []; hemiID = [];
        for c = 1:numel(ctKeys)
            ctname = ctKeys{c};
            S_all  = hemiDataByType.(ctname);
            if isempty(S_all), continue; end

            S = S_all(strcmp({S_all.sex}, sexKey));
            if isempty(S), continue; end

            densVec = vertcat(S.density);           % [nHemi × nSub]
            nHemi   = size(densVec,1);

            % Select the column corresponding to this subregion
            y_c = densVec(:, rIdx);                 % density for this subregion
            hemi_c = strings(nHemi,1);
            for h = 1:nHemi
                hemi_c(h) = sprintf('%s_%s', S(h).brainID, S(h).hemi);
            end

            y      = [y; y_c];
            hemiID = [hemiID; hemi_c];
            ct     = [ct; repmat({ctname}, nHemi, 1)];
        end

        % Clean up
        ok = isfinite(y);
        y = y(ok);
        if isempty(y)
            fprintf('Subregion %s (sex %s): no hemispheres; skipping.\n', region, sexKey);
            continue;
        end

        ct     = categorical(ct(ok));
        hemiID = categorical(hemiID(ok));

        % Mixed-effects model
        T = table(y, ct, hemiID, 'VariableNames', {'Density','CellType','HemiID'});
        mdl = fitlme(T, 'Density ~ 1 + CellType + (1|HemiID)', 'DummyVarCoding','reference');
        A   = anova(mdl);
        disp(A)

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
        allRegionP(rIdx, sSex) = p_sub;
        fprintf('Subregion %s (sex %s): omnibus CellType effect p = %.3g\n', region, sexKey, p_sub);

        % --- Post-hoc pairwise between cell types (per subregion, this sex) ---
        ctCats = categories(ct);
        K      = numel(ctCats);

        % Fixed effects and covariance (coerce to numeric if needed)
        beta = fixedEffects(mdl);                    % p×1
        CovB = mdl.CoefficientCovariance;
        if ~isnumeric(CovB)
            if istable(CovB)
                CovB = table2array(CovB);
            elseif isa(CovB,'dataset') || isa(CovB,'titleddataset')
                CovB = table2array(dataset2table(CovB));
            else
                CovB = double(CovB);   % last resort
            end
        end

        coefNames = string(mdl.CoefficientNames);    % e.g., "(Intercept)","CellType_TH",...

        % Map each CellType level to its coefficient index (empty -> reference level)
        lvl2idx = containers.Map('KeyType','char','ValueType','any');
        for k = 1:K
            lvl = ctCats{k};
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

        fprintf('[%s, sex %s] LME pairwise (CellType) via fixed-effect contrasts (Sidak):\n', region, sexKey);
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
end

% Cross-region FDR, per sex
for sSex = 1:numel(sexKeys)
    sexKey = sexKeys{sSex};
    pcol = allRegionP(:, sSex);
    ok = isfinite(pcol);
    if ~any(ok), continue; end

    pvals = pcol(ok);
    [~, qvals] = fdr_bh_local(pvals);
    fprintf('\n== CellType effect across subregions (BH-FDR) — sex %s ==\n', sexKey);
    for k = 1:numel(pvals)
        fprintf('  %s: raw p = %.3g, q = %.3g\n', allRegionNames{k}, pvals(k), qvals(k));
    end
end


%% 6. Quick visualization to check if clusterVol is 1 hemisphere or whole brain
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

%% 7. QC: plot one AP plane of detections colored by 5/7/12/15

ctQC    = 'PV';        % choose CT to inspect (e.g. 'TH','SST','PV')
brainQC = '10222';     % choose brain ID to inspect

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
                [N,~,bin] = histcounts(ap, edges); %#ok<ASGLU>
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


%% helper functions
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