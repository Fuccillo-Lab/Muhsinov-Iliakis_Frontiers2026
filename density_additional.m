%% DENSITY PROCESSING AND PLOT
%This code was used to generate supplemental figures in Muhsinov et al.
%(2026). 
%Input data files available on Zenodo: https://zenodo.org/records/18685812

%% 0. Setup and Parameters
addpath('C:\Users\walki\Box\grp-psom-fuccillo-lab\Interneuron Anatomy\Scripts');
workingDir = 'C:\Users\walki\Box\grp-psom-fuccillo-lab\Interneuron Anatomy\Brain Analysis\Brains';
cd(workingDir);

voxelSize = 150; %can change (um/pixel)
voxelVol_mm3 = (voxelSize / 1000)^3;   % (0.15 mm)^3 = 0.003375 mm^3

Bregma = 5400; % Approximation of Bregma Z coord (posterior from anterior)
filtSigma = 0.5; % Gaussian smoothing sigma (in *voxel* units of the grid)
Nboot = 1000;
sliceSpacing = voxelSize/1000;   % mm per slice
% Midline in CCFv3 (um). In Allen CCFv3, ML spans 0 to 11400 with midline ~5700 (MBF classifies switch at 5695).

midline = 5700;
width   = 2*midline;   % 11400

% Block 2 — Multi–cell-type ingest (TH/SST/PV), hemisphere split, and right-flip

REGION = 'Caudoputamen';

% Columns in the classified sheets that hold ML/DV/AP in µm
COL_ML = 8;
COL_DV = 9;
COL_AP = 10;

csvFiles = dir(fullfile(workingDir, '*_4clusters.xlsx'));  % all brains/all types

% Only these types are processed
ctKeys = {'TH','SST','PV'};   % add 'ChAT' etc when ready
sexKeys = {'F','M'};

%% 1. Import Data
% Containers
hemiDataByType   = struct();
allDataByType    = struct();
fileIndexByType  = struct();
for k = 1:numel(ctKeys)
    hemiDataByType.(ctKeys{k})  = struct('file',{},'brainID',{},'sex',{},'hemi',{},'coords',{});
    allDataByType.(ctKeys{k})   = [];
    fileIndexByType.(ctKeys{k}) = struct('file',{},'brainID',{},'sex',{});
end

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


for f = csvFiles'
    [~, base, ~] = fileparts(f.name);

    % Try raw pattern, then classified pattern
    m = regexp(base, pat_raw, 'names', 'once');
    if isempty(m)
        m = regexp(base, pat_classified, 'names', 'once');
    end

    if isempty(m)
        warning('PARSE FAIL (skipping): %s', f.name);
        continue
    end

    ct      = normCellType(m.ct);       % e.g. PV/SST/TH; assumes this exists
    brainID = m.id;
    sex     = upper(m.sex);

    if isempty(ct) || ~isfield(allDataByType, ct)
        warning('Unsupported/Excluded cell type "%s" (skipping file): %s', m.ct, f.name);
        continue
    end

    T = readtable(fullfile(f.folder, f.name));

    % Region filter
    if ~ismember('Region', T.Properties.VariableNames)
        warning('%s: no "Region" column → skipping', f.name);
        continue
    end
    T = T(strcmpi(T.Region, REGION), :);
    if isempty(T)
        continue
    end

    % insideMask = T.inside_mask;
    % if ~islogical(insideMask)
    %     insideMask = logical(insideMask);   % handles 0/1 or similar
    % end
    %
    % % Keep only rows in REGION AND inside the mask
    % keep = strcmpi(T.Region, REGION) & insideMask;
    % T    = T(keep, :);

    % Hemisphere split (unflipped)
    TL = T(T{:,COL_ML} <  midline, :);   % left
    TR = T(T{:,COL_ML} >= midline, :);   % right

    if ~isempty(TL)
        hemiDataByType.(ct)(end+1) = struct( ...
            'file',    fullfile(f.folder, f.name), ...
            'brainID', brainID, ...
            'sex',     sex, ...
            'hemi',    'L', ...
            'coords',  double(TL{:, [COL_ML COL_DV COL_AP]}) );
    end
    if ~isempty(TR)
        hemiDataByType.(ct)(end+1) = struct( ...
            'file',    fullfile(f.folder, f.name), ...
            'brainID', brainID, ...
            'sex',     sex, ...
            'hemi',    'R', ...
            'coords',  double(TR{:, [COL_ML COL_DV COL_AP]}) );
    end

    % Right-hemisphere aggregate (flip left to right)
    % Requires "midline" and "width" already defined in your Block 1.
    ml  = double(T{:,COL_ML});
    idx = (ml < midline);
    ml(idx) = -ml(idx) + width;   % reflect left points across midline
    coordsR = [ml, double(T{:,COL_DV}), double(T{:,COL_AP})];

    allDataByType.(ct) = [allDataByType.(ct); coordsR]; %vertconcatenating
    fileIndexByType.(ct)(end+1) = struct( ...
        'file',    fullfile(f.folder, f.name), ...
        'brainID', brainID, ...
        'sex',     sex );
end

% Combined across included types (right-flipped)
allData = [];
for k = 1:numel(ctKeys)
    allData = [allData; allDataByType.(ctKeys{k})]; %#ok<AGROW>
end

% Summary
for k = 1:numel(ctKeys)
    nFiles = numel(fileIndexByType.(ctKeys{k}));
    nPts   = size(allDataByType.(ctKeys{k}),1);
    if nFiles > 0
        fprintf('%s: %d file(s), %d point(s)\n', ctKeys{k}, nFiles, nPts);
    end
end
M = sum(cellfun(@(ct) numel(hemiDataByType.(ct)), ctKeys));
fprintf('Built %d hemisphere buckets across included types.\n', M);

%% 2a. Build voxels
% Convert to voxel‐indices and define integer bin edges
% Per-type voxel indices
voxelIdxByType = struct();
for k = 1:numel(ctKeys)
    ct = ctKeys{k};
    coords_um = allDataByType.(ct);
    if isempty(coords_um)
        voxelIdxByType.(ct) = zeros(0,3);
    else
        voxelIdxByType.(ct) = round(coords_um / voxelSize);  % [iX iY iZ]
    end
end

% Global grid extents using all types combined (keeps heatmaps aligned)
voxelIdxAll = round(allData / voxelSize);   % each row now [iX iY iZ] in voxel units

imin = min(voxelIdxAll, [], 1);   % [iXmin, iYmin, iZmin]
imax = max(voxelIdxAll, [], 1);   % [iXmax, iYmax, iZmax]

% Bin edges *between* integer voxel centers (bin width = 1 index = 1 voxel)
pts = {
    (imin(1)-0.5) : 1 : (imax(1)+0.5),    % x-edges
    (imin(2)-0.5) : 1 : (imax(2)+0.5),    % y-edges
    (imin(3)-0.5) : 1 : (imax(3)+0.5)     % z-edges
    };

gridSz = [numel(pts{1})-1, numel(pts{2})-1, numel(pts{3})-1]; % We have N edges in pts, so N-1 bins

% Per-type, per-sex accumulators (for bootstrap running sums/means/vars)
accumByType = struct();
for k = 1:numel(ctKeys)
    ct = ctKeys{k};
    for s = 1:numel(sexKeys)
        sex = sexKeys{s};  % 'F' or 'M'
        accumByType.(ct).(sex) = struct( ...
            'runningSum',          zeros(gridSz), ...
            'runningSum_unsmooth', zeros(gridSz));
    end
end

% (Optional) quick summary
fprintf('Grid size: [%d %d %d]\n', gridSz(1), gridSz(2), gridSz(3));
for k = 1:numel(ctKeys)
    ct = ctKeys{k};
    fprintf('%s: %d points, voxelIdx size %dx3\n', ct, size(allDataByType.(ct),1), size(voxelIdxByType.(ct),1));
end

%% 2b. Exclude voxels outside caudoputamen
% Build a striatum mask in CCFv3 (25 µm) and map to your 150 µm grid
% Start with the local 150-µm grid (gridSz) and shift it by imin to get **global** 150-µm indices,
% then convert those indices to physical positions in microns by multiplying by voxelSize.
% Map these um coordinates to the **25-µm** atlas by taking round(…/25) (nearest neighbor), and clamp the resulting indices to the atlas bounds.
% Next, use `ndgrid`/`sub2ind` to resample the 25-µm striatum mask (`M25`) onto your 150-µm grid—this yields `mask150`, which flags whether each analysis voxel lies inside striatum.
% Optionally, apply a **core-only** step (L∞ erosion) so that a full 150-um cube centered at any kept voxel is guaranteed to be entirely within tissue.
%M25 = the actual 25 µm striatum mask (3-D logical). ml_25/dv_25/ap_25 = index lookups into M25. mask150 = the mask on your 150 µm grid after resampling M25.

cpFile = 'structure_672.nrrd';  % 0 = outside, 1 = inside
% --- Load NRRD and make a logical mask at 25 µm ---
V = nrrdread(cpFile);          % V: [z y x] = [AP DV ML]

% Reorder to (ML, DV, AP) to match the rest of your code
M25 = permute(V > 0, [3 2 1]); % now size(M25) = [ML, DV, AP]

% (Optional) grab dims if you want them explicitly
ML = size(M25,1);
DV = size(M25,2);
AP = size(M25,3);

% --- Build bin-center coordinates (um) for the shared 150 µm grid ---
% Use the global grid defined earlier (from all types combined)
nX = gridSz(1); nY = gridSz(2); nZ = gridSz(3);

ml_idx150 = imin(1) : imax(1);
dv_idx150 = imin(2) : imax(2);
ap_idx150 = imin(3) : imax(3);

ml_um = ml_idx150 * voxelSize;
dv_um = dv_idx150 * voxelSize;
ap_um = ap_idx150 * voxelSize;

% If centers align between grids (both CCFv3-style), use ROUND:
ml_25 = round(ml_um / 25) + 1;
dv_25 = round(dv_um / 25) + 1;
ap_25 = round(ap_um / 25) + 1;

% Guard against out-of-bounds
%ml_25 = min(max(ml_25, 1), size(M25,1));
%dv_25 = min(max(dv_25, 1), size(M25,2));
%ap_25 = min(max(ap_25, 1), size(M25,3));

assert(all(ml_25>=1 & ml_25<=size(M25,1)));
assert(all(dv_25>=1 & dv_25<=size(M25,2)));
assert(all(ap_25>=1 & ap_25<=size(M25,3)));

% --- Sample the 25 µm mask onto the 150 µm grid (nearest-neighbor) ---
[IX, IY, IZ] = ndgrid(ml_25, dv_25, ap_25);  % order matches M25: (ML, DV, AP)
mask150 = M25(sub2ind(size(M25), IX, IY, IZ));  % [nX nY nZ] logical



% Core-only mask (optional)
rim_um   = 50;                 % 50 µm "inset" from boundary
rim_vox  = ceil(rim_um / 25);  % = 2 voxels at 25 µm

d_inf = bwdist(~M25, 'chessboard');   % L∞ distance in 25-µm voxels
core25 = d_inf >= rim_vox;            % at least rim_um from boundary

mask150_core = core25(sub2ind(size(M25), IX, IY, IZ));   % [nX nY nZ] logical

% Diagnostics and final assign
coverage_core = nnz(mask150_core) / max(nnz(mask150), 1);
fprintf('Core coverage (>= %d vox = %dum): %.1f%%\n', rim_vox, rim_vox*25, 100*coverage_core);

mask150_precore = mask150;
mask150 = mask150_core;   % comment this out if you want the full CP instead of core-only

% Ventricle buffer: exclude voxels < 175 µm (Chebyshev) from ventricles

ventFile = 'Ventricle_Mask.nrrd';   % 25 µm ventricle mask (0/1)
Vv = nrrdread(ventFile);            % [z y x] = [AP DV ML]
vent25 = permute(Vv > 0, [3 2 1]);  % [ML DV AP], logical

% Chebyshev (L∞) distance in 25-µm voxels to nearest ventricle voxel
dvent_inf = bwdist(vent25, 'chessboard');   % 0 inside ventricle

ventRim_um  = 175;                           % desired ventricle buffer (µm)
ventRim_vox = ceil(ventRim_um / 25);

% Keep only voxels at least 175 µm away from ventricles
ventFar25 = dvent_inf >= ventRim_vox;       % 1 = ≥ 175 µm from ventricle

% Resample ventricle buffer to 150-µm grid usig same IX,IY,IZ map
ventFar150 = ventFar25(sub2ind(size(ventFar25), IX, IY, IZ));  % [nX nY nZ] logical

% Combine CP core mask with ventricle buffer
mask150 = mask150 & ventFar150; % can remove this by commenting out

%% 3. Map cells to voxels: Bootstrap → Bin → Smooth (per cell type, per sex)

% per-type, per-sex outputs
px_bootByType = struct();   % px_bootByType.(ct).F, px_bootByType.(ct).M  [nX x Nboot]
py_bootByType = struct();   % py_bootByType.(ct).F, ...
pz_bootByType = struct();   % pz_bootByType.(ct).F, ...

meanVol_unsmoothByType = struct();  % meanVol_unsmoothByType.(ct).F / .M
meanVolByType          = struct();  % meanVolByType.(ct).F / .M

for c = 1:numel(ctKeys)
    ct     = ctKeys{c};
    H_all  = hemiDataByType.(ct);   % all hemispheres for this cell type

    % Initialize per-sex containers for this CT
    px_bootByType.(ct)          = struct();
    py_bootByType.(ct)          = struct();
    pz_bootByType.(ct)          = struct();
    meanVol_unsmoothByType.(ct) = struct();
    meanVolByType.(ct)          = struct();

    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};   % 'F' or 'M'

        % Subset hemisphere buckets for this sex
        if isempty(H_all)
            H = H_all;
        else
            H = H_all(strcmp({H_all.sex}, sexKey));
        end
        M = numel(H);  % number of hemispheres for this ct & sex

        if M == 0
            % No hemispheres of this sex for this CT → store empties and skip
            fprintf('[%s] No hemispheres for sex %s; skipping bootstrap.\n', ct, sexKey);
            px_bootByType.(ct).(sexKey)          = [];
            py_bootByType.(ct).(sexKey)          = [];
            pz_bootByType.(ct).(sexKey)          = [];
            meanVol_unsmoothByType.(ct).(sexKey) = [];
            meanVolByType.(ct).(sexKey)          = [];
            continue;
        end

        % Guard: accumulators for this ct/sex should be zeroed
        if any(accumByType.(ct).(sexKey).runningSum(:)) || ...
                any(accumByType.(ct).(sexKey).runningSum_unsmooth(:))
            error('accumByType.%s.%s accumulators are non-zero. Clear or re-create them before bootstrapping.', ...
                ct, sexKey);
        end

        % Allocate per-sex bootstrap profiles for this CT
        px_bootByType.(ct).(sexKey) = zeros(nX, Nboot);
        py_bootByType.(ct).(sexKey) = zeros(nY, Nboot);
        pz_bootByType.(ct).(sexKey) = zeros(nZ, Nboot);

        for b = 1:Nboot
            % --- Resample hemisphere buckets with replacement (within sex) ---
            if M == 0
                samplePts = [];   % already handled above, but keep guard
            else
                pick = randsample(M, M, true);  % bootstrap hemispheres
                samplePts = zeros(0,3);
                for kH = 1:M
                    C = H(pick(kH)).coords;     % [ML DV AP] (µm), unflipped
                    if isempty(C), continue; end
                    if H(pick(kH)).hemi == 'L'  % flip left → right
                        C(:,1) = -C(:,1) + width;  % ML' = -ML + 2*midline
                    end
                    samplePts = [samplePts; C];
                end
            end

            % Bin to 3D histogram on the shared grid
            if isempty(samplePts)
                counts = zeros(gridSz);
            else
                sampIdx = round(samplePts / voxelSize);    % → voxel indices [iX iY iZ]
                counts  = histcnd(sampIdx, pts);           % 1-voxel bins
            end

            % Harmonize histcnd output: verify extra edge planes are empty before trimming
            edgesSz = [numel(pts{1}), numel(pts{2}), numel(pts{3})];

            if isequal(size(counts), edgesSz)
                % sums over the 3 "last-edge" planes
                sx = sum(sum(sum(counts(end,   :,   :))));
                sy = sum(sum(sum(counts(:,   end,   :))));
                sz = sum(sum(sum(counts(:,     :, end))));
                % Histcnd returns an array with size equal to number of edges,
                % hence we must drop the empty plane
                if sx==0 && sy==0 && sz==0
                    % Safe to drop the empty overflow planes → convert edges→bins
                    counts = counts(1:end-1, 1:end-1, 1:end-1);
                else
                    error('Non-zero last-edge planes [sx=%d sy=%d sz=%d]...', sx, sy, sz);
                end
            end

            % Final safety assert: counts now matches your gridSz
            assert(isequal(size(counts), gridSz), ...
                'counts size %s != gridSz %s', mat2str(size(counts)), mat2str(gridSz));

            % Calculate DENSITY, NOT PROBABILITY
            tot = sum(counts(:));
            if tot == 0
                pVol_unsmooth = zeros(gridSz);
            else
                % Divide by number of hemispheres of this sex (M) to get per-hemisphere density
                pVol_unsmooth = (counts ./ M) ./ voxelVol_mm3;
            end
            %
            % % --- 1D profiles from UNSMOOTHED density volume (mask outside as NaN) ---
            % P_b = pVol_unsmooth;
            % P_b(~mask150) = NaN;
            %
            % px_bootByType.(ct).(sexKey)(:, b) = squeeze( mean(P_b, [2 3], 'omitnan') );  % ML
            % py_bootByType.(ct).(sexKey)(:, b) = squeeze( mean(P_b, [1 3], 'omitnan') );  % DV
            % pz_bootByType.(ct).(sexKey)(:, b) = squeeze( mean(P_b, [1 2], 'omitnan') );  % AP

            % Update UNSMOOTHED accumulators for this ct/sex
            accumByType.(ct).(sexKey).runningSum_unsmooth = ...
                accumByType.(ct).(sexKey).runningSum_unsmooth + pVol_unsmooth;

            % SMOOTHED density volume (Gaussian sigma in voxel units)
            if tot == 0
                pVol_smooth = zeros(gridSz);
            else
                pVol_smooth = imgaussfilt3(pVol_unsmooth, filtSigma);
            end
            %  1D profiles from SMOOTHEDdensity volume (mask outside as NaN)
            P_b = pVol_smooth;
            P_b(~mask150) = NaN;
            %
            px_bootByType.(ct).(sexKey)(:, b) = squeeze( mean(P_b, [2 3], 'omitnan') );  % ML
            py_bootByType.(ct).(sexKey)(:, b) = squeeze( mean(P_b, [1 3], 'omitnan') );  % DV
            pz_bootByType.(ct).(sexKey)(:, b) = squeeze( mean(P_b, [1 2], 'omitnan') );  % AP

            % --- Update SMOOTHED accumulators for this ct/sex ---
            accumByType.(ct).(sexKey).runningSum = ...
                accumByType.(ct).(sexKey).runningSum + pVol_smooth;
        end  % bootstrap loop

        % --- Per-sex mean volumes (unsmoothed + smoothed) for this cell type ---
        meanVol_unsmoothByType.(ct).(sexKey) = ...
            accumByType.(ct).(sexKey).runningSum_unsmooth / Nboot;

        meanVolByType.(ct).(sexKey) = ...
            accumByType.(ct).(sexKey).runningSum / Nboot;

        % (Optionally) print a little summary
        fprintf('[%s, sex=%s] Bootstraps complete: %d hemis, %d boots.\n', ...
            ct, sexKey, M, Nboot);
    end  % sex loop
end  % ct loop


%% 4. Mask per-type volumes and print pre/post sums (per sex)

for c = 1:numel(ctKeys)
    ct = ctKeys{c};

    if ~isfield(meanVolByType, ct)
        continue;
    end

    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};

        if ~isfield(meanVolByType.(ct), sexKey) || ...
                isempty(meanVolByType.(ct).(sexKey))
            continue;
        end

        mv_uns = meanVol_unsmoothByType.(ct).(sexKey);
        mv_smo = meanVolByType.(ct).(sexKey);

        % Pre-mask sums
        s_uns_pre = sum(mv_uns(:)) * voxelVol_mm3;
        s_smo_pre = sum(mv_smo(:)) * voxelVol_mm3;
        fprintf('[%s %s] Pre-mask:  sum(meanVol_unsmooth)=%.6f,  sum(meanVol)=%.6f\n', ...
            ct, sexKey, s_uns_pre, s_smo_pre);

        % Apply mask (use NaN so masked voxels don’t affect means/plots with omitnan)
        mv_uns(~mask150) = NaN;
        mv_smo(~mask150) = NaN;

        % Write back
        meanVol_unsmoothByType.(ct).(sexKey) = mv_uns;
        meanVolByType.(ct).(sexKey)          = mv_smo;

        % Post-mask sums (inside striatum)
        s_uns_post = nansum(mv_uns(:)) * voxelVol_mm3;
        s_smo_post = nansum(mv_smo(:)) * voxelVol_mm3;
        fprintf('[%s %s] Post-mask: nansum(meanVol_unsmooth)=%.6f, nansum(meanVol)=%.6f\n', ...
            ct, sexKey, s_uns_post, s_smo_post);
    end
end

% ========================================================================
% ADDED BLOCK: Bootstrap -> Bin -> Smooth (per cell type, per hemisphere)
% Collapsed across sex to preserve power.
% Left hemispheres are mirrored into right-hemisphere space before binning,
% so L and R can be directly compared in a common coordinate frame.
% ========================================================================

% Per-type, per-hemisphere outputs
px_bootByType_hemi = struct();   % px_bootByType_hemi.(ct).L / .R  [nX x Nboot]
py_bootByType_hemi = struct();   % py_bootByType_hemi.(ct).L / .R  [nY x Nboot]
pz_bootByType_hemi = struct();   % pz_bootByType_hemi.(ct).L / .R  [nZ x Nboot]

meanVol_unsmoothByType_hemi = struct();  % meanVol_unsmoothByType_hemi.(ct).L / .R
meanVolByType_hemi          = struct();  % meanVolByType_hemi.(ct).L / .R

% Separate accumulators for hemisphere analysis
accumByType_hemi = struct();
for k = 1:numel(ctKeys)
    ct = ctKeys{k};
    for h = 1:numel(hemiKeys)
        hemi = hemiKeys{h};  % 'L' or 'R'
        accumByType_hemi.(ct).(hemi) = struct( ...
            'runningSum',          zeros(gridSz), ...
            'runningSum_unsmooth', zeros(gridSz));
    end
end

%% 5. Bootstrap in a separate group by hemi - sex-agnostic; to allow for comparison of densities across different hemispheres in statistical analyses
hemiKeys = {'L','R'};
% Initialize hemisphere accumulators
accumByType_hemi = struct();
% Pre-initialize hemisphere accumulators for every ct x hemi combination
for k = 1:numel(ctKeys)
    ct = ctKeys{k};
    for h = 1:numel(hemiKeys)
        hemi = hemiKeys{h};
        accumByType_hemi.(ct).(hemi) = struct( ...
            'runningSum', zeros(gridSz), ...
            'runningSum_unsmooth', zeros(gridSz) );
    end
end
for c = 1:numel(ctKeys)
    ct    = ctKeys{c};
    H_all = hemiDataByType.(ct);   % all hemisphere buckets for this cell type

    % Initialize containers for this cell type
    px_bootByType_hemi.(ct)          = struct();
    py_bootByType_hemi.(ct)          = struct();
    pz_bootByType_hemi.(ct)          = struct();
    meanVol_unsmoothByType_hemi.(ct) = struct();
    meanVolByType_hemi.(ct)          = struct();

    for h = 1:numel(hemiKeys)
        hemiKey = hemiKeys{h};   % 'L' or 'R'

        % Subset by hemisphere ONLY (collapsed across sex)
        if isempty(H_all)
            H = H_all;
        else
            H = H_all(strcmp({H_all.hemi}, hemiKey));
        end
        M = numel(H);  % number of hemisphere buckets for this ct & hemi

        if M == 0
            fprintf('[HEMI %s] No hemisphere buckets for %s; skipping bootstrap.\n', ct, hemiKey);
            px_bootByType_hemi.(ct).(hemiKey)          = [];
            py_bootByType_hemi.(ct).(hemiKey)          = [];
            pz_bootByType_hemi.(ct).(hemiKey)          = [];
            meanVol_unsmoothByType_hemi.(ct).(hemiKey) = [];
            meanVolByType_hemi.(ct).(hemiKey)          = [];
            continue;
        end

        % Safety check
        if any(accumByType_hemi.(ct).(hemiKey).runningSum(:)) || ...
           any(accumByType_hemi.(ct).(hemiKey).runningSum_unsmooth(:))
            error('accumByType_hemi.%s.%s accumulators are non-zero. Clear or re-create them before bootstrapping.', ...
                ct, hemiKey);
        end

        % Allocate bootstrap profile storage
        px_bootByType_hemi.(ct).(hemiKey) = zeros(nX, Nboot);
        py_bootByType_hemi.(ct).(hemiKey) = zeros(nY, Nboot);
        pz_bootByType_hemi.(ct).(hemiKey) = zeros(nZ, Nboot);

        for b = 1:Nboot
            % --- Resample hemisphere buckets with replacement ---
            pick = randsample(M, M, true);
            samplePts = zeros(0,3);

            for kH = 1:M
                C = H(pick(kH)).coords;   % [ML DV AP] in um
                if isempty(C), continue; end

                % Mirror LEFT hemisphere into RIGHT-space for direct comparison
                if hemiKey == 'L'
                    C(:,1) = -C(:,1) + width;
                end

                samplePts = [samplePts; C]; %#ok<AGROW>
            end

            % Bin to 3D histogram on shared grid
            if isempty(samplePts)
                counts = zeros(gridSz);
            else
                sampIdx = round(samplePts / voxelSize);
                counts  = histcnd(sampIdx, pts);
            end

            % Harmonize histcnd output
            edgesSz = [numel(pts{1}), numel(pts{2}), numel(pts{3})];

            if isequal(size(counts), edgesSz)
                sx = sum(sum(sum(counts(end, :, :))));
                sy = sum(sum(sum(counts(:, end, :))));
                sz = sum(sum(sum(counts(:, :, end))));
                if sx==0 && sy==0 && sz==0
                    counts = counts(1:end-1, 1:end-1, 1:end-1);
                else
                    error('Non-zero last-edge planes [sx=%d sy=%d sz=%d] in hemisphere bootstrap.', sx, sy, sz);
                end
            end

            assert(isequal(size(counts), gridSz), ...
                'counts size %s != gridSz %s', mat2str(size(counts)), mat2str(gridSz));

            % Per-hemisphere density
            tot = sum(counts(:));
            if tot == 0
                pVol_unsmooth = zeros(gridSz);
            else
                pVol_unsmooth = (counts ./ M) ./ voxelVol_mm3;
            end

            % Update unsmoothed running sum
            accumByType_hemi.(ct).(hemiKey).runningSum_unsmooth = ...
                accumByType_hemi.(ct).(hemiKey).runningSum_unsmooth + pVol_unsmooth;

            % Smooth
            if tot == 0
                pVol_smooth = zeros(gridSz);
            else
                pVol_smooth = imgaussfilt3(pVol_unsmooth, filtSigma);
            end

            % 1D profiles from smoothed density volume, masked outside striatum
            P_b = pVol_smooth;
            P_b(~mask150) = NaN;

            px_bootByType_hemi.(ct).(hemiKey)(:, b) = squeeze(mean(P_b, [2 3], 'omitnan')); % ML
            py_bootByType_hemi.(ct).(hemiKey)(:, b) = squeeze(mean(P_b, [1 3], 'omitnan')); % DV
            pz_bootByType_hemi.(ct).(hemiKey)(:, b) = squeeze(mean(P_b, [1 2], 'omitnan')); % AP

            % Update smoothed running sum
            accumByType_hemi.(ct).(hemiKey).runningSum = ...
                accumByType_hemi.(ct).(hemiKey).runningSum + pVol_smooth;
        end

        % Mean volumes across bootstrap iterations
        meanVol_unsmoothByType_hemi.(ct).(hemiKey) = ...
            accumByType_hemi.(ct).(hemiKey).runningSum_unsmooth / Nboot;

        meanVolByType_hemi.(ct).(hemiKey) = ...
            accumByType_hemi.(ct).(hemiKey).runningSum / Nboot;

        fprintf('[HEMI %s, %s] Bootstraps complete: %d hemisphere buckets, %d boots.\n', ...
            ct, hemiKey, M, Nboot);
    end
end


%% 6. Sanity checks & quick visuals (per cell type)

% Pick AP slice to inspect: use BregmaIdx if available; else nearest to Bregma (µm)
if exist('BregmaIdx','var')
    sliceIdx = BregmaIdx;
else
    % ap_idx150 are the global 150-µm voxel indices along AP (built earlier in mask block)
    targetAPidx = round(Bregma / voxelSize);
    [~, sliceIdx] = min(abs(ap_idx150 - targetAPidx));
end

% 0) One-time assertion: mask must match grid & any one meanVol (any ct, any sex)
refCT  = '';
refSex = '';

for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    if ~isfield(meanVolByType, ct), continue; end
    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};
        if isfield(meanVolByType.(ct), sexKey) && ...
                ~isempty(meanVolByType.(ct).(sexKey))
            refCT  = ct;
            refSex = sexKey;
            break;
        end
    end
    if ~isempty(refCT)
        break;
    end
end

assert(~isempty(refCT), 'No meanVol found in meanVolByType.');
assert(isequal(size(mask150), size(meanVolByType.(refCT).(refSex))), ...
    'mask150 and meanVol size mismatch');

% 1) View the mask alone (same for all CTs)
figure('Name','Striatum mask (mask150) — coronal slice');
imagesc(squeeze(mask150(:,:,sliceIdx))');   % LM × DV view
set(gca,'XDir','reverse'); axis image off;
colormap(gray);
title(sprintf('mask150 at AP slice %d', sliceIdx));

% 2) Overlay mask outline on each CT × sex’s mean volume for alignment check
for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    if ~isfield(meanVolByType, ct), continue; end

    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};
        if ~isfield(meanVolByType.(ct), sexKey) || ...
                isempty(meanVolByType.(ct).(sexKey))
            continue;
        end

        MV = meanVolByType.(ct).(sexKey);

        figure('Name', sprintf('[%s %s] MeanVol with striatum outline', ct, sexKey));
        img = squeeze(MV(:,:,sliceIdx))';
        imagesc(img, 'AlphaData', ~isnan(img));    % hide NaNs outside mask if present
        set(gca,'XDir','reverse'); axis image off; colormap(flipud(gray(256)));
        hold on;
        edge = bwperim(squeeze(mask150(:,:,sliceIdx))');
        [hY,hX] = find(edge);
        plot(hX, hY, '.', 'MarkerSize', 4);
        title(sprintf('[%s %s] AP slice %d with striatum outline', ...
            ct, sexKey, sliceIdx));
    end
end

% 3) Quick coverage stats for this slice (mask-only; same for all CTs)
insideFrac = nnz(squeeze(mask150(:,:,sliceIdx))) / numel(squeeze(mask150(:,:,sliceIdx)));
fprintf('Slice %d: %0.1f%% of bins are inside the striatum mask.\n', sliceIdx, insideFrac*100);

%% 7. 1D Density profiles across anatomical axes (per cell type, per sex) with shaded SEM

% Axes in mm (shared for all CTs; built from previously computed bin centers)
x_mm = (ml_um - midline) / 1000;      % ML: mm from midline
y_mm = (dv_um - dv_um(1)) / 1000;     % DV: mm from dorsal-most plane
z_mm = -(ap_um - Bregma) / 1000;      % AP: mm from Bregma (anterior positive)

% Outputs: nested by ct and sex
densityXByType = struct();  semXByType = struct();
densityYByType = struct();  semYByType = struct();
densityZByType = struct();  semZByType = struct();

% Colors per sex (keep consistent across plots)
sexColor.F = [0.8500 0.3250 0.0980];  % female: orange-ish
sexColor.M = [0.0000 0.4470 0.7410];  % male: blue
lightFactor = 0.5;                    % for SEM shading

for c = 1:numel(ctKeys)
    ct = ctKeys{c};

    % --- Compute densities/SEMs per sex for this CT ---
    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};

        if ~isfield(px_bootByType, ct) || ...
                ~isfield(px_bootByType.(ct), sexKey) || ...
                isempty(px_bootByType.(ct).(sexKey))
            continue;
        end

        pxB = px_bootByType.(ct).(sexKey);  % [nX x Nboot]
        pyB = py_bootByType.(ct).(sexKey);  % [nY x Nboot]
        pzB = pz_bootByType.(ct).(sexKey);  % [nZ x Nboot]

        densityXByType.(ct).(sexKey) = nanmean(pxB, 2);
        densityYByType.(ct).(sexKey) = nanmean(pyB, 2);
        densityZByType.(ct).(sexKey) = nanmean(pzB, 2);

        semXByType.(ct).(sexKey) = nanstd(pxB, [], 2);
        semYByType.(ct).(sexKey) = nanstd(pyB, [], 2);
        semZByType.(ct).(sexKey) = nanstd(pzB, [], 2);
    end

    % If no sex had data, skip CT
    if ~isfield(densityXByType, ct)
        fprintf('[%s] No bootstrap profiles for any sex; skipping density plots.\n', ct);
        continue;
    end

    % ========= ML axis =========
    figML = figure('Name', sprintf('[%s] Density — ML (by sex)', ct));
    hold on; box on; grid on;

    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};
        if ~isfield(densityXByType.(ct), sexKey) || ...
                isempty(densityXByType.(ct).(sexKey))
            continue;
        end

        baseColor = sexColor.(sexKey);
        semColor  = baseColor + (1 - baseColor)*lightFactor;

        x_full  = x_mm(:);
        y_full  = densityXByType.(ct).(sexKey)(:);
        se_full = semXByType.(ct).(sexKey)(:);

        assert(numel(x_full) == numel(y_full) && numel(y_full) == numel(se_full), ...
            '[%s ML %s] Length mismatch between x, mean, and SEM.', ct, sexKey);

        valid = isfinite(y_full) & isfinite(se_full);
        x = x_full(valid);
        y = y_full(valid);
        se = se_full(valid);
        if isempty(x), continue; end

        upper = y + se;

        % Fill under mean
        xx = [x; flipud(x)];
        yy = [y; zeros(size(y))];
        fill(xx, yy, baseColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        % SEM band (mean → mean+SE)
        xx2 = [x; flipud(x)];
        yy2 = [upper; flipud(y)];
        fill(xx2, yy2, semColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        % Lines
        plot(x, y,     'Color', baseColor, 'LineWidth', 2);
        plot(x, upper, 'Color', semColor,  'LineWidth', 1.2);
    end

    ylim([0, 1750]);
    xlabel('ML (mm from midline)');
    ylabel('Detection density');
    title(sprintf('[%s] ML density with SEM band (by sex)', ct));
    legend(sexKeys, 'Location','northeast');
    hold off;

    % ========= DV axis =========
    figDV = figure('Name', sprintf('[%s] Density — DV (by sex)', ct));
    hold on; box on; grid on;

    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};
        if ~isfield(densityYByType.(ct), sexKey) || ...
                isempty(densityYByType.(ct).(sexKey))
            continue;
        end

        baseColor = sexColor.(sexKey);
        semColor  = baseColor + (1 - baseColor)*lightFactor;

        x_full  = y_mm(:);
        y_full  = densityYByType.(ct).(sexKey)(:);
        se_full = semYByType.(ct).(sexKey)(:);

        assert(numel(x_full) == numel(y_full) && numel(y_full) == numel(se_full), ...
            '[%s DV %s] Length mismatch between x, mean, and SEM.', ct, sexKey);

        valid = isfinite(y_full) & isfinite(se_full);
        x = x_full(valid);
        y = y_full(valid);
        se = se_full(valid);
        if isempty(x), continue; end

        upper = y + se;

        xx = [x; flipud(x)];
        yy = [y; zeros(size(y))];
        fill(xx, yy, baseColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        xx2 = [x; flipud(x)];
        yy2 = [upper; flipud(y)];
        fill(xx2, yy2, semColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        plot(x, y,     'Color', baseColor, 'LineWidth', 2);
        plot(x, upper, 'Color', semColor,  'LineWidth', 1.2);
    end

    ylim([0, 2500]);
    xlabel('DV (mm from dorsal)');
    ylabel('Detection density');
    title(sprintf('[%s] DV density with SEM band (by sex)', ct));
    legend(sexKeys, 'Location','northeast');
    hold off;

    % ========= AP axis =========
    figAP = figure('Name', sprintf('[%s] Density — AP (by sex)', ct));
    hold on; box on; grid on;

    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};
        if ~isfield(densityZByType.(ct), sexKey) || ...
                isempty(densityZByType.(ct).(sexKey))
            continue;
        end

        baseColor = sexColor.(sexKey);
        semColor  = baseColor + (1 - baseColor)*lightFactor;

        x_full  = z_mm(:);
        y_full  = densityZByType.(ct).(sexKey)(:);
        se_full = semZByType.(ct).(sexKey)(:);

        assert(numel(x_full) == numel(y_full) && numel(y_full) == numel(se_full), ...
            '[%s AP %s] Length mismatch between x, mean, and SEM.', ct, sexKey);

        valid = isfinite(y_full) & isfinite(se_full);
        x = x_full(valid);
        y = y_full(valid);
        se = se_full(valid);
        if isempty(x), continue; end

        upper = y + se;

        xx = [x; flipud(x)];
        yy = [y; zeros(size(y))];
        fill(xx, yy, baseColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        xx2 = [x; flipud(x)];
        yy2 = [upper; flipud(y)];
        fill(xx2, yy2, semColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        plot(x, y,     'Color', baseColor, 'LineWidth', 2);
        plot(x, upper, 'Color', semColor,  'LineWidth', 1.2);
    end

    ylim([0, 2250]);
    set(gca, 'XDir','reverse');  % preserve AP convention
    xlabel('AP (mm from Bregma)');
    ylabel('Detection density');
    title(sprintf('[%s] AP density with SEM band (by sex)', ct));
    legend(sexKeys, 'Location','northeast');
    hold off;
end

%% 8. 1D Density plots per sex - tiled layout

% Axes in mm (shared for all CTs; built from previously computed bin centers)
x_mm = (ml_um - midline) / 1000;      % ML: mm from midline
y_mm = (dv_um - dv_um(1)) / 1000;     % DV: mm from dorsal-most plane
z_mm = -(ap_um - Bregma) / 1000;      % AP: mm from Bregma (anterior positive)

% Outputs: nested by ct and sex
densityXByType = struct();  semXByType = struct();
densityYByType = struct();  semYByType = struct();
densityZByType = struct();  semZByType = struct();

% Colors per sex (keep consistent across plots)
sexColor.F = [0.8500 0.3250 0.0980];  % female: orange-ish
sexColor.M = [0.0000 0.4470 0.7410];  % male: blue
lightFactor = 0.5;                    % for SEM shading
ctKeys = {"SST", "PV", "TH"};
legendPrinted = 0;
figure; tiledlayout(3,3);

for c = 1:numel(ctKeys)
    ct = ctKeys{c};

    % --- Compute densities/SEMs per sex for this CT ---
    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};

        if ~isfield(px_bootByType, ct) || ...
                ~isfield(px_bootByType.(ct), sexKey) || ...
                isempty(px_bootByType.(ct).(sexKey))
            continue;
        end

        pxB = px_bootByType.(ct).(sexKey);  % [nX x Nboot]
        pyB = py_bootByType.(ct).(sexKey);  % [nY x Nboot]
        pzB = pz_bootByType.(ct).(sexKey);  % [nZ x Nboot]

        densityXByType.(ct).(sexKey) = nanmean(pxB, 2);
        densityYByType.(ct).(sexKey) = nanmean(pyB, 2);
        densityZByType.(ct).(sexKey) = nanmean(pzB, 2);

        semXByType.(ct).(sexKey) = nanstd(pxB, [], 2);
        semYByType.(ct).(sexKey) = nanstd(pyB, [], 2);
        semZByType.(ct).(sexKey) = nanstd(pzB, [], 2);
    end

    % If no sex had data, skip CT
    if ~isfield(densityXByType, ct)
        fprintf('[%s] No bootstrap profiles for any sex; skipping density plots.\n', ct);
        continue;
    end

    % ========= AP axis =========
    nexttile;
    title('[%s] Density — AP (by sex)', ct);
    hold on; box on; grid on;

    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};
        if ~isfield(densityZByType.(ct), sexKey) || ...
                isempty(densityZByType.(ct).(sexKey))
            continue;
        end

        baseColor = sexColor.(sexKey);
        semColor  = baseColor + (1 - baseColor)*lightFactor;

        x_full  = z_mm(:);
        y_full  = densityZByType.(ct).(sexKey)(:);
        se_full = semZByType.(ct).(sexKey)(:);

        assert(numel(x_full) == numel(y_full) && numel(y_full) == numel(se_full), ...
            '[%s AP %s] Length mismatch between x, mean, and SEM.', ct, sexKey);

        valid = isfinite(y_full) & isfinite(se_full);
        x = x_full(valid);
        y = y_full(valid);
        se = se_full(valid);
        if isempty(x), continue; end

        upper = y + se;

        xx = [x; flipud(x)];
        yy = [y; zeros(size(y))];
        fill(xx, yy, baseColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        xx2 = [x; flipud(x)];
        yy2 = [upper; flipud(y)];
        fill(xx2, yy2, semColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        plot(x, y,     'Color', baseColor, 'LineWidth', 2);
        plot(x, upper, 'Color', semColor,  'LineWidth', 1.2);
    end

    ylim([0, 2250]);
    set(gca, 'XDir','reverse');  % preserve AP convention
    xlabel('AP (mm from Bregma)');
    ylabel('Detection density');
    title(sprintf('[%s] AP density with SE band (by sex)', ct));
    hold off;

    % ========= ML axis =========
    nexttile;
    title('[%s] Density — ML (by sex)', ct);
    hold on; box on; grid on;

    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};
        if ~isfield(densityXByType.(ct), sexKey) || ...
                isempty(densityXByType.(ct).(sexKey))
            continue;
        end

        baseColor = sexColor.(sexKey);
        semColor  = baseColor + (1 - baseColor)*lightFactor;

        x_full  = x_mm(:);
        y_full  = densityXByType.(ct).(sexKey)(:);
        se_full = semXByType.(ct).(sexKey)(:);

        assert(numel(x_full) == numel(y_full) && numel(y_full) == numel(se_full), ...
            '[%s ML %s] Length mismatch between x, mean, and SEM.', ct, sexKey);

        valid = isfinite(y_full) & isfinite(se_full);
        x = x_full(valid);
        y = y_full(valid);
        se = se_full(valid);
        if isempty(x), continue; end

        upper = y + se;

        % Fill under mean
        xx = [x; flipud(x)];
        yy = [y; zeros(size(y))];
        fill(xx, yy, baseColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        % SEM band (mean → mean+SE)
        xx2 = [x; flipud(x)];
        yy2 = [upper; flipud(y)];
        fill(xx2, yy2, semColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        % Lines
        plot(x, y,     'Color', baseColor, 'LineWidth', 2);
        plot(x, upper, 'Color', semColor,  'LineWidth', 1.2);
    end

    ylim([0, 1750]);
    xlabel('ML (mm from midline)');
    ylabel('Detection density');
    title(sprintf('[%s] ML density with SE band (by sex)', ct));
    hold off;

    % ========= DV axis =========
    nexttile;
    title('[%s] Density — DV (by sex)', ct);
    hold on; box on; grid on;

    for s = 1:numel(sexKeys)
        sexKey = sexKeys{s};
        if ~isfield(densityYByType.(ct), sexKey) || ...
                isempty(densityYByType.(ct).(sexKey))
            continue;
        end

        baseColor = sexColor.(sexKey);
        semColor  = baseColor + (1 - baseColor)*lightFactor;

        x_full  = y_mm(:);
        y_full  = densityYByType.(ct).(sexKey)(:);
        se_full = semYByType.(ct).(sexKey)(:);

        assert(numel(x_full) == numel(y_full) && numel(y_full) == numel(se_full), ...
            '[%s DV %s] Length mismatch between x, mean, and SEM.', ct, sexKey);

        valid = isfinite(y_full) & isfinite(se_full);
        x = x_full(valid);
        y = y_full(valid);
        se = se_full(valid);
        if isempty(x), continue; end

        upper = y + se;

        xx = [x; flipud(x)];
        yy = [y; zeros(size(y))];
        fill(xx, yy, baseColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        xx2 = [x; flipud(x)];
        yy2 = [upper; flipud(y)];
        fill(xx2, yy2, semColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        plot(x, y,     'Color', baseColor, 'LineWidth', 2);
        plot(x, upper, 'Color', semColor,  'LineWidth', 1.2);
    end

    ylim([0, 2500]);
    xlabel('DV (mm from dorsal)');
    ylabel('Detection density');
    title(sprintf('[%s] DV density with SE band (by sex)', ct));
    if legendPrinted == 0
        legendPrinted = 1;
        text(0.98, 0.95, 'Female', ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'Color', [0.85 0.325 0.098], ...
            'FontWeight', 'bold');

        text(0.98, 0.83, 'Male', ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'Color', [0 0.447 0.741], ...
            'FontWeight', 'bold');
    end
    hold off;


end




%% 9. Mixed-effects regression per cell type and axis — BY HEMISPHERE (1D density profiles)

% Model: Y ~ X + (1|hemiID)       (random intercept per hemisphere)
% Option: Y ~ X + (1+X|hemiID)    (random intercept + slope; may be singular with few hemis)
%
% This block:
%   - builds per-hemisphere 1D density profiles (ML/DV/AP),
%   - stacks them into long tables,
%   - fits LME per axis for each cell type,
%   - plots data + fixed-effect line (+ optional spaghetti for random slopes).
%
% Requires:
%   ctKeys, hemiDataByType, voxelSize, voxelVol_mm3, pts, gridSz, mask150,
%   ml_um, dv_um, ap_um, Bregma already defined.

useRandomSlope = false;
showHemiSlopes = false;

lmeModels = struct();     % optional
ctKeys = {"SST", "PV", "TH"};
colF = [0.85 0.325 0.098];
colM = [0 0.447 0.741];
legendPrinted = 0;

figure; tiledlayout(3,3);
for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    H  = hemiDataByType.(ct);            % hemisphere buckets for this cell type

    if isempty(H)
        fprintf('[%s] No hemisphere buckets; skipping LME.\n', ct);
        continue;
    end

    % Build per-hemisphere density distributions
    byHemi = struct('hemiID', {}, 'sex', {}, 'px', {}, 'py', {}, 'pz', {});
    hCount = 0;

    for k = 1:numel(H)
        C = H(k).coords;                    % [ML DV AP] in um, unflipped
        if isempty(C), continue; end

        C2 = C;
        if H(k).hemi == 'L'
            % reflect left → right into common CCF space
            C2(:,1) = -C2(:,1) + width;     % ML' = -ML + 2*midline
        end
        % C2 is all coordinates
        % Have to use this helper function to build density distributions PER HEMISPHERE
        [px, py, pz] = profile_one_hemi(C2, voxelSize, pts, gridSz, mask150, voxelVol_mm3);
        if isempty(px) || isempty(py) || isempty(pz)
            continue;
        end

        hCount = hCount + 1;
        byHemi(hCount).hemiID = sprintf('%s_%s', H(k).brainID, H(k).hemi);  % e.g. "04_R"
        byHemi(hCount).sex    = H(k).sex;   % 'F' or 'M'
        byHemi(hCount).px     = px;   % column vector, ML density profile
        byHemi(hCount).py     = py;   % DV density profile
        byHemi(hCount).pz     = pz;   % AP density profile
    end

    if isempty(byHemi)
        fprintf('[%s] No per-hemisphere profiles after filtering; skipping LME.\n', ct);
        continue;
    end

    A  = struct2table(byHemi);
    nH = height(A);

    % ---------- Axis vectors (in mm) ----------
    % Recreate from earlier bin-center coords to avoid accidental edits.
    x_mm_axis = (ml_um - midline) / 1000;        % ML: mm from midline
    y_mm_axis = (dv_um - dv_um(1)) / 1000;       % DV: mm from dorsal-most plane
    z_mm_axis = -(ap_um - Bregma) / 1000;        % AP: mm from Bregma (ant +)

    % Long-format tables: one row = one bin of one hemisphere

    % ML
    Yx = []; Xx = []; IDx = strings(0,1); Sexx = strings(0,1);
    for i = 1:nH
        vx = getcol(A.px, i); vx = vx(:);
        assert(numel(vx) == numel(x_mm_axis), ...
            '[%s ML] Length mismatch between hemisphere %d profile and x_mm_axis.', ct, i);
    
        Yx   = [Yx; vx];
        Xx   = [Xx; x_mm_axis(:)];
        IDx  = [IDx; repmat(string(A.hemiID{i}), numel(vx), 1)];
        Sexx = [Sexx; repmat(string(A.sex{i}),    numel(vx), 1)];
    end
    
    Tx = table(Yx, Xx, categorical(IDx), categorical(Sexx, {'F','M'}), ...
        'VariableNames', {'Y','X','hemiID','sex'});


    % DV
   Yy = []; Xy = []; IDy = strings(0,1); Sexy = strings(0,1);
    for i = 1:nH
        vy = getcol(A.py, i); vy = vy(:);
        assert(numel(vy) == numel(y_mm_axis), ...
            '[%s DV] Length mismatch between hemisphere %d profile and y_mm_axis.', ct, i);
    
        Yy   = [Yy; vy];
        Xy   = [Xy; y_mm_axis(:)];
        IDy  = [IDy; repmat(string(A.hemiID{i}), numel(vy), 1)];
        Sexy = [Sexy; repmat(string(A.sex{i}),    numel(vy), 1)];
    end
    
    Ty = table(Yy, Xy, categorical(IDy), categorical(Sexy, {'F','M'}), ...
        'VariableNames', {'Y','X','hemiID','sex'});


    % AP
    Yz = []; Xz = []; IDz = strings(0,1); Sexz = strings(0,1);
    for i = 1:nH
        vz = getcol(A.pz, i); vz = vz(:);
        assert(numel(vz) == numel(z_mm_axis), ...
            '[%s AP] Length mismatch between hemisphere %d profile and z_mm_axis.', ct, i);
    
        Yz   = [Yz; vz];
        Xz   = [Xz; z_mm_axis(:)];
        IDz  = [IDz; repmat(string(A.hemiID{i}), numel(vz), 1)];
        Sexz = [Sexz; repmat(string(A.sex{i}),    numel(vz), 1)];
    end
    
    Tz = table(Yz, Xz, categorical(IDz), categorical(Sexz, {'F','M'}), ...
        'VariableNames', {'Y','X','hemiID','sex'});


    % Keep raw X for plotting
    Xx_raw = Xx;
    Xy_raw = Xy;
    Xz_raw = Xz;

    % ---------- Fit LME models ----------
    if useRandomSlope
        formula = 'Y ~ X + (1 + X|hemiID)';
        modelStr = 'Y ~ X + (1+X|hemiID)';
    else
        formula  = 'Y ~ X*sex + (1|hemiID)';
        modelStr = 'Y ~ X*sex + (1|hemiID)';


    end

    % Random slope fit
        mdlX = fitlme(Tx, formula);
        mdlY = fitlme(Ty, formula);
        mdlZ = fitlme(Tz, formula);


    % ---------- Report fixed-effect slope β_X ----------
    CX = mdlX.Coefficients; rowX = strcmp(CX.Name,'X');
    CY = mdlY.Coefficients; rowY = strcmp(CY.Name,'X');
    CZ = mdlZ.Coefficients; rowZ = strcmp(CZ.Name,'X');

    fprintf('\n[%s] Mixed-effects by hemisphere: %s\n', ct, modelStr);
    fprintf(' ML: beta_X = %.4g (SE=%.4g), p=%.3g, nHemi=%d, nRows=%d\n', ...
        CX.Estimate(rowX), CX.SE(rowX), CX.pValue(rowX), nH, height(Tx));
    fprintf(' DV: beta_X = %.4g (SE=%.4g), p=%.3g, nHemi=%d, nRows=%d\n', ...
        CY.Estimate(rowY), CY.SE(rowY), CY.pValue(rowY), nH, height(Ty));
    fprintf(' AP: beta_X = %.4g (SE=%.4g), p=%.3g, nHemi=%d, nRows=%d\n', ...
        CZ.Estimate(rowZ), CZ.SE(rowZ), CZ.pValue(rowZ), nH, height(Tz));

    % Store models
    lmeModels.(ct).ML = mdlX;
    lmeModels.(ct).DV = mdlY;
    lmeModels.(ct).AP = mdlZ;

    % ---------- Plots: data + fixed-effect line + optional spaghetti ----------

    annFmt = '\\beta = %.3g (SE=%.3g), p = %.3g';

    % --- AP ---
    nexttile; hold on; box on; grid on;

    isF = (Tz.sex == 'F');
    isM = (Tz.sex == 'M');
    
    plot(Tz.X(isF), Tz.Y(isF), '.', 'MarkerSize', 6, 'Color', colF);
    plot(Tz.X(isM), Tz.Y(isM), '.', 'MarkerSize', 6, 'Color', colM);
    
    zg = z_mm_axis(:);
    
    % fixed-effect (population) predictions by sex
    [yhatF, yhatM] = predictBySexFixed(mdlZ, Tz, zg);
    
    plot(zg, yhatF, 'LineWidth', 2, 'Color', colF);
    plot(zg, yhatM, 'LineWidth', 2, 'Color', colM);
    
    set(gca,'XDir','reverse');
    xlabel('AP (mm from Bregma)');
    ylabel('Density (cells / mm^3 / hemisphere)');
    title(sprintf('[%s] AP: %s', ct, modelStr));
    


    % --- ML ---
    nexttile; hold on; box on; grid on;

    isF = (Tx.sex == 'F');
    isM = (Tx.sex == 'M');
    
    plot(Tx.X(isF), Tx.Y(isF), '.', 'MarkerSize', 6, 'Color', colF);
    plot(Tx.X(isM), Tx.Y(isM), '.', 'MarkerSize', 6, 'Color', colM);
    
    xg = x_mm_axis(:);
    
    [yhatF, yhatM] = predictBySexFixed(mdlX, Tx, xg);
    
    plot(xg, yhatF, 'LineWidth', 2, 'Color', colF);
    plot(xg, yhatM, 'LineWidth', 2, 'Color', colM);
    
    xlabel('ML (mm from midline)');
    ylabel('Density (cells / mm^3 / hemisphere)');
    title(sprintf('[%s] ML: %s', ct, modelStr));


    % --- DV ---
    nexttile; hold on; box on; grid on;

    isF = (Ty.sex == 'F');
    isM = (Ty.sex == 'M');
    
    plot(Ty.X(isF), Ty.Y(isF), '.', 'MarkerSize', 6, 'Color', colF);
    plot(Ty.X(isM), Ty.Y(isM), '.', 'MarkerSize', 6, 'Color', colM);
    
    yg = y_mm_axis(:);
    
    [yhatF, yhatM] = predictBySexFixed(mdlY, Ty, yg);
    
    plot(yg, yhatF, 'LineWidth', 2, 'Color', colF);
    plot(yg, yhatM, 'LineWidth', 2, 'Color', colM);
    
    xlabel('DV (mm from dorsal)');
    ylabel('Density (cells / mm^3 / hemisphere)');
    title(sprintf('[%s] DV: %s', ct, modelStr));
    legendPrinted = legendPrinted + 1;
    if legendPrinted == 2
        legendPrinted = 3;
        text(0.98, 0.95, 'Female', 'Units','normalized', 'HorizontalAlignment','right', ...
             'VerticalAlignment','top', 'Color', colF, 'FontWeight','bold');
        text(0.98, 0.83, 'Male',   'Units','normalized', 'HorizontalAlignment','right', ...
             'VerticalAlignment','top', 'Color', colM, 'FontWeight','bold');
    end



    
end

%% 10. RIGHT VS. LEFT DENSITY COMPARISON (Mixed-effects regression per cell type and axis; 1D density profiles)

% Model: Y ~ X*hemi + (1|brainID)
%
% This block:
%   - builds per-hemisphere 1D density profiles (ML/DV/AP),
%   - mirrors LEFT hemispheres into RIGHT-space before profiling,
%   - stacks them into long tables,
%   - fits LME per axis for each cell type,
%   - plots data + fixed-effect lines for L and R.
%
% Requires:
%   ctKeys, hemiDataByType, voxelSize, voxelVol_mm3, pts, gridSz, mask150,
%   ml_um, dv_um, ap_um, Bregma, midline, width already defined.

useRandomSlope = false;

lmeModels = struct();
ctKeys = {"SST", "PV", "TH"};

colL = [0.85 0.325 0.098];
colR = [0 0.447 0.741];
legendPrinted = 0;

figure; tiledlayout(3,3);

for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    H  = hemiDataByType.(ct);

    if isempty(H)
        fprintf('[%s] No hemisphere buckets; skipping LME.\n', ct);
        continue;
    end

    % Build per-hemisphere density distributions
    byHemi = struct('brainID', {}, 'hemiID', {}, 'hemi', {}, 'px', {}, 'py', {}, 'pz', {});
    hCount = 0;

    for k = 1:numel(H)
        C = H(k).coords;   % [ML DV AP] in um, unflipped
        if isempty(C), continue; end

        C2 = C;
        if H(k).hemi == 'L'
            % reflect left -> right into common space
            C2(:,1) = -C2(:,1) + width;
        end

        [px, py, pz] = profile_one_hemi(C2, voxelSize, pts, gridSz, mask150, voxelVol_mm3);
        if isempty(px) || isempty(py) || isempty(pz)
            continue;
        end

        hCount = hCount + 1;
        byHemi(hCount).brainID = string(H(k).brainID);                      % e.g. "02"
        byHemi(hCount).hemiID  = string(sprintf('%s_%s', H(k).brainID, H(k).hemi)); % e.g. "02_L"
        byHemi(hCount).hemi    = string(H(k).hemi);                         % 'L' or 'R'
        byHemi(hCount).px      = px;
        byHemi(hCount).py      = py;
        byHemi(hCount).pz      = pz;
    end

    if isempty(byHemi)
        fprintf('[%s] No per-hemisphere profiles after filtering; skipping LME.\n', ct);
        continue;
    end

    A  = struct2table(byHemi);
    nH = height(A);
    nBrains = numel(unique(A.brainID));

    % ---------- Axis vectors (in mm) ----------
    x_mm_axis = (ml_um - midline) / 1000;   % ML: mm from midline
    y_mm_axis = (dv_um - dv_um(1)) / 1000;  % DV: mm from dorsal-most plane
    z_mm_axis = -(ap_um - Bregma) / 1000;   % AP: mm from Bregma (anterior positive)

    % ---------- Long-format tables ----------

    % ML
    Yx = []; Xx = []; BrainX = strings(0,1); HemiX = strings(0,1); HemiIDX = strings(0,1);
    for i = 1:nH
        vx = getcol(A.px, i); vx = vx(:);
        assert(numel(vx) == numel(x_mm_axis), ...
            '[%s ML] Length mismatch between hemisphere %d profile and x_mm_axis.', ct, i);

        Yx      = [Yx; vx];
        Xx      = [Xx; x_mm_axis(:)];
        BrainX  = [BrainX; repmat(string(A.brainID(i)), numel(vx), 1)];
        HemiX   = [HemiX;  repmat(string(A.hemi(i)),    numel(vx), 1)];
        HemiIDX = [HemiIDX; repmat(string(A.hemiID(i)), numel(vx), 1)];
    end

    Tx = table(Yx, Xx, categorical(BrainX), categorical(HemiX, {'L','R'}), categorical(HemiIDX), ...
        'VariableNames', {'Y','X','brainID','hemi','hemiID'});

    % DV
    Yy = []; Xy = []; BrainY = strings(0,1); HemiY = strings(0,1); HemiIDY = strings(0,1);
    for i = 1:nH
        vy = getcol(A.py, i); vy = vy(:);
        assert(numel(vy) == numel(y_mm_axis), ...
            '[%s DV] Length mismatch between hemisphere %d profile and y_mm_axis.', ct, i);

        Yy      = [Yy; vy];
        Xy      = [Xy; y_mm_axis(:)];
        BrainY  = [BrainY; repmat(string(A.brainID(i)), numel(vy), 1)];
        HemiY   = [HemiY;  repmat(string(A.hemi(i)),    numel(vy), 1)];
        HemiIDY = [HemiIDY; repmat(string(A.hemiID(i)), numel(vy), 1)];
    end

    Ty = table(Yy, Xy, categorical(BrainY), categorical(HemiY, {'L','R'}), categorical(HemiIDY), ...
        'VariableNames', {'Y','X','brainID','hemi','hemiID'});

    % AP
    Yz = []; Xz = []; BrainZ = strings(0,1); HemiZ = strings(0,1); HemiIDZ = strings(0,1);
    for i = 1:nH
        vz = getcol(A.pz, i); vz = vz(:);
        assert(numel(vz) == numel(z_mm_axis), ...
            '[%s AP] Length mismatch between hemisphere %d profile and z_mm_axis.', ct, i);

        Yz      = [Yz; vz];
        Xz      = [Xz; z_mm_axis(:)];
        BrainZ  = [BrainZ; repmat(string(A.brainID(i)), numel(vz), 1)];
        HemiZ   = [HemiZ;  repmat(string(A.hemi(i)),    numel(vz), 1)];
        HemiIDZ = [HemiIDZ; repmat(string(A.hemiID(i)), numel(vz), 1)];
    end

    Tz = table(Yz, Xz, categorical(BrainZ), categorical(HemiZ, {'L','R'}), categorical(HemiIDZ), ...
        'VariableNames', {'Y','X','brainID','hemi','hemiID'});

    % ---------- Fit LME models ----------
    if useRandomSlope
        formula  = 'Y ~ X*hemi + (1 + X|brainID)';
        modelStr = 'Y ~ X*hemi + (1+X|brainID)';
    else
        formula  = 'Y ~ X*hemi + (1|brainID)';
        modelStr = 'Y ~ X*hemi + (1|brainID)';
    end

    mdlX = fitlme(Tx, formula);
    mdlY = fitlme(Ty, formula);
    mdlZ = fitlme(Tz, formula);

    % ---------- Report exact fixed effects ----------
    CX = mdlX.Coefficients;
    CY = mdlY.Coefficients;
    CZ = mdlZ.Coefficients;

    rowX_X      = strcmp(CX.Name, 'X');
    rowX_hemiR  = strcmp(CX.Name, 'hemi_R');
    rowX_XhemiR = strcmp(CX.Name, 'X:hemi_R');

    rowY_X      = strcmp(CY.Name, 'X');
    rowY_hemiR  = strcmp(CY.Name, 'hemi_R');
    rowY_XhemiR = strcmp(CY.Name, 'X:hemi_R');

    rowZ_X      = strcmp(CZ.Name, 'X');
    rowZ_hemiR  = strcmp(CZ.Name, 'hemi_R');
    rowZ_XhemiR = strcmp(CZ.Name, 'X:hemi_R');

    fprintf('\n[%s] Mixed-effects by hemisphere: %s\n', ct, modelStr);
    fprintf(' nHemi=%d, nBrains=%d\n', nH, nBrains);

    fprintf(' ML: X = %.4g (SE=%.4g), p=%.3g; hemi_R = %.4g (SE=%.4g), p=%.3g; X:hemi_R = %.4g (SE=%.4g), p=%.3g; nRows=%d\n', ...
        CX.Estimate(rowX_X),      CX.SE(rowX_X),      CX.pValue(rowX_X), ...
        CX.Estimate(rowX_hemiR),  CX.SE(rowX_hemiR),  CX.pValue(rowX_hemiR), ...
        CX.Estimate(rowX_XhemiR), CX.SE(rowX_XhemiR), CX.pValue(rowX_XhemiR), ...
        height(Tx));

    fprintf(' DV: X = %.4g (SE=%.4g), p=%.3g; hemi_R = %.4g (SE=%.4g), p=%.3g; X:hemi_R = %.4g (SE=%.4g), p=%.3g; nRows=%d\n', ...
        CY.Estimate(rowY_X),      CY.SE(rowY_X),      CY.pValue(rowY_X), ...
        CY.Estimate(rowY_hemiR),  CY.SE(rowY_hemiR),  CY.pValue(rowY_hemiR), ...
        CY.Estimate(rowY_XhemiR), CY.SE(rowY_XhemiR), CY.pValue(rowY_XhemiR), ...
        height(Ty));

    fprintf(' AP: X = %.4g (SE=%.4g), p=%.3g; hemi_R = %.4g (SE=%.4g), p=%.3g; X:hemi_R = %.4g (SE=%.4g), p=%.3g; nRows=%d\n', ...
        CZ.Estimate(rowZ_X),      CZ.SE(rowZ_X),      CZ.pValue(rowZ_X), ...
        CZ.Estimate(rowZ_hemiR),  CZ.SE(rowZ_hemiR),  CZ.pValue(rowZ_hemiR), ...
        CZ.Estimate(rowZ_XhemiR), CZ.SE(rowZ_XhemiR), CZ.pValue(rowZ_XhemiR), ...
        height(Tz));

    % Store models
    lmeModels.(ct).ML = mdlX;
    lmeModels.(ct).DV = mdlY;
    lmeModels.(ct).AP = mdlZ;

    % ---------- Plots: AP ----------
    nexttile; hold on; box on; grid on;

    isL = (Tz.hemi == 'L');
    isR = (Tz.hemi == 'R');

    plot(Tz.X(isL), Tz.Y(isL), '.', 'MarkerSize', 6, 'Color', colL);
    plot(Tz.X(isR), Tz.Y(isR), '.', 'MarkerSize', 6, 'Color', colR);

    zg = z_mm_axis(:);
    [yhatL, yhatR] = predictByHemiFixed(mdlZ, Tz, zg);

    plot(zg, yhatL, 'LineWidth', 2, 'Color', colL);
    plot(zg, yhatR, 'LineWidth', 2, 'Color', colR);

    set(gca,'XDir','reverse');
    xlabel('AP (mm from Bregma)');
    ylabel('Density (cells / mm^3 / hemisphere)');
    title(sprintf('[%s] AP: %s', ct, modelStr));

    % ---------- Plots: ML ----------
    nexttile; hold on; box on; grid on;

    isL = (Tx.hemi == 'L');
    isR = (Tx.hemi == 'R');

    plot(Tx.X(isL), Tx.Y(isL), '.', 'MarkerSize', 6, 'Color', colL);
    plot(Tx.X(isR), Tx.Y(isR), '.', 'MarkerSize', 6, 'Color', colR);

    xg = x_mm_axis(:);
    [yhatL, yhatR] = predictByHemiFixed(mdlX, Tx, xg);

    plot(xg, yhatL, 'LineWidth', 2, 'Color', colL);
    plot(xg, yhatR, 'LineWidth', 2, 'Color', colR);

    xlabel('ML (mm from midline)');
    ylabel('Density (cells / mm^3 / hemisphere)');
    title(sprintf('[%s] ML: %s', ct, modelStr));

    % ---------- Plots: DV ----------
    nexttile; hold on; box on; grid on;

    isL = (Ty.hemi == 'L');
    isR = (Ty.hemi == 'R');

    plot(Ty.X(isL), Ty.Y(isL), '.', 'MarkerSize', 6, 'Color', colL);
    plot(Ty.X(isR), Ty.Y(isR), '.', 'MarkerSize', 6, 'Color', colR);

    yg = y_mm_axis(:);
    [yhatL, yhatR] = predictByHemiFixed(mdlY, Ty, yg);

    plot(yg, yhatL, 'LineWidth', 2, 'Color', colL);
    plot(yg, yhatR, 'LineWidth', 2, 'Color', colR);

    xlabel('DV (mm from dorsal)');
    ylabel('Density (cells / mm^3 / hemisphere)');
    title(sprintf('[%s] DV: %s', ct, modelStr));

    legendPrinted = legendPrinted + 1;
    if legendPrinted == 2
        legendPrinted = 3;
        text(0.98, 0.95, 'Left', 'Units','normalized', 'HorizontalAlignment','right', ...
             'VerticalAlignment','top', 'Color', colL, 'FontWeight','bold');
        text(0.98, 0.83, 'Right', 'Units','normalized', 'HorizontalAlignment','right', ...
             'VerticalAlignment','top', 'Color', colR, 'FontWeight','bold');
    end
end

%% Helper functions

function [px, py, pz] = profile_one_hemi(C, voxelSize, pts, gridSz, mask150, voxelVol_mm3)
% Bin a hemisphere’s detections to the shared grid and return 1D DENSITY profiles.
% C: [ML DV AP] in µm, already flipped into right hemisphere if needed.

if isempty(C)
    px = []; py = []; pz = [];
    return;
end

% Map to voxel indices, same as histcnd logic above
sampIdx = round(C / voxelSize);
ix = discretize(sampIdx(:,1), pts{1});
iy = discretize(sampIdx(:,2), pts{2});
iz = discretize(sampIdx(:,3), pts{3});

keep = isfinite(ix) & isfinite(iy) & isfinite(iz);
if ~any(keep)
    px = []; py = []; pz = [];
    return;
end

subs   = [ix(keep), iy(keep), iz(keep)];
counts = accumarray(subs, 1, gridSz);    % raw counts in each voxel

if ~any(counts(:))
    px = []; py = []; pz = [];
    return;
end

% Convert to density: cells / mm^3 / hemisphere
% (single hemisphere’s data here → no division by M)
P = counts ./ voxelVol_mm3;

% Mask outside CP
P(~mask150) = NaN;

% 1D marginal density profiles (mean across other dims, ignoring NaNs)
px = squeeze( mean(P, [2 3], 'omitnan') ); px = px(:);
py = squeeze( mean(P, [1 3], 'omitnan') ); py = py(:);
pz = squeeze( mean(P, [1 2], 'omitnan') ); pz = pz(:);
end

function v = getcol(M, i)
% Fetch the i-th profile from table column M.
% Handles both cell arrays and numeric matrices.
if iscell(M)
    v = M{i};
else
    % Assume rows = hemisphere, cols = bins
    v = M(i, :);
end
end

function plot_spaghetti_lines(mdl, tbl_axis, x_grid, color_rgb)
hemi_ids = categories(tbl_axis.hemiID);
if isempty(hemi_ids)
    warning('plot_spaghetti_lines: no hemisphere categories; skipping.');
    return;
end

n_pts = numel(x_grid);

for h = 1:numel(hemi_ids)
    this_id = hemi_ids{h};

    hemi_col = categorical(repmat({this_id}, n_pts, 1), hemi_ids);
    pred_tbl = table(x_grid, hemi_col, 'VariableNames', {'X','hemiID'});

    y_hat = predict(mdl, pred_tbl, 'Conditional', true);
    plot(x_grid, y_hat, '-', 'LineWidth', 0.8, 'Color', color_rgb);
end
end

function [yhatF, yhatM] = predictBySexFixed(mdl, T, x_grid)
% Fixed-effect predictions for sex = F and sex = M on x_grid.
% Uses a dummy hemiID level (required by fitlme) but Conditional=false
% ensures population-level (fixed-effects) prediction.

n = numel(x_grid);

hemiCats = categories(T.hemiID);
dummyHemi = categorical(repmat(hemiCats(1), n, 1), hemiCats);

sexCats = categories(T.sex);
predF = table(x_grid(:), dummyHemi, categorical(repmat("F", n, 1), sexCats), ...
    'VariableNames', {'X','hemiID','sex'});
predM = table(x_grid(:), dummyHemi, categorical(repmat("M", n, 1), sexCats), ...
    'VariableNames', {'X','hemiID','sex'});

yhatF = predict(mdl, predF, 'Conditional', false);
yhatM = predict(mdl, predM, 'Conditional', false);
end

function [yhatL, yhatR] = predictByHemiFixed(mdl, T, x_grid)
% Fixed-effect predictions for hemi = L and hemi = R on x_grid.
% Uses a dummy brainID level from the existing categories; Conditional=false
% returns population-level (fixed-effects) prediction.

n = numel(x_grid);

brainCats = categories(T.brainID);
dummyBrain = categorical(repmat(brainCats(1), n, 1), brainCats);

hemiCats = categories(T.hemi);

predL = table(x_grid(:), dummyBrain, categorical(repmat("L", n, 1), hemiCats), ...
    'VariableNames', {'X','brainID','hemi'});

predR = table(x_grid(:), dummyBrain, categorical(repmat("R", n, 1), hemiCats), ...
    'VariableNames', {'X','brainID','hemi'});

yhatL = predict(mdl, predL, 'Conditional', false);
yhatR = predict(mdl, predR, 'Conditional', false);
end





% ---- helper ----
function s = normCellType(raw)
r = upper(regexprep(raw,'\+',''));  % drop '+' if present
switch r
    case 'CHAT', s = 'ChAT';
    case {'TH','SST','PV'}, s = r;
    otherwise, s = '';  % excluded/unsupported
end
end

% ---- Local helper handles both scalar and arrays of Bar objects ----
function setBarProps(bh)
if numel(bh) > 1
    set(bh, 'EdgeColor','none');   % older MATLAB: set across array
    set(bh, 'FaceAlpha',0.9);
else
    bh.EdgeColor = 'none';
    bh.FaceAlpha = 0.9;
end
end
