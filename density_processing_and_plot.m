%% DENSITY PROCESSING AND PLOT
%This code was used to generate Figures 2 and 3 of Muhsinov et al. (2026).
%Input data files available on Zenodo: https://zenodo.org/records/18685812

%% 0. Setup and Parameters
addpath('C:\Users\walki\Box\grp-psom-fuccillo-lab\Interneuron Anatomy\Scripts');
workingDir = 'C:\Users\walki\Box\grp-psom-fuccillo-lab\Interneuron Anatomy\Brain Analysis\Brains';
cd(workingDir);
voxelSize = 150; % (um/pixel)
voxelVol_mm3 = (voxelSize / 1000)^3;   % (0.15 mm)^3 = 0.003375 mm^3

Bregma = 5400; % Approximation of Bregma Z coord (posterior from anterior)
filtSigma = 0.5; % Gaussian smoothing sigma (in voxel units of the grid)
Nboot = 1000;
sliceSpacing = voxelSize/1000;   % mm per slice
% Midline in CCFv3 (um). In Allen CCFv3, ML spans 0 to 11400 with midline ~5700

midline = 5700;
width   = 2*midline;   % 11400

% Block 2: Multi–cell-type ingest (TH/SST/PV), hemisphere split, and right-flip

REGION = 'Caudoputamen';

% Columns in the classified sheets that hold ML/DV/AP in microns
COL_ML = 8;
COL_DV = 9;
COL_AP = 10;

csvFiles = dir(fullfile(workingDir, '*_4clusters.xlsx'));  % all brains/all types

ctKeys = {'TH','SST','PV'};   % cell types

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

% Filename patterns for raw MBF output csv files:
% TH_04_F_points or TH_04_F_new_points
pat_raw = '^(?<ct>[A-Za-z\+]+)_(?<id>[0-9A-Za-z]+)_(?<sex>[MFmf])(?:_new)?_points$';

% Filename patterns for subregion classified xlsx files:
% PV_10222_F_points_classified_4clusters or  TH_F02_new_points_classified_4clusters
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

    ct      = normCellType(m.ct);       % PV/SST/TH; assumes this exists
    brainID = m.id;
    sex     = upper(m.sex);

    if isempty(ct) || ~isfield(allDataByType, ct)
        warning('Unsupported/Excluded cell type "%s" (skipping file): %s', m.ct, f.name);
        continue
    end

    T = readtable(fullfile(f.folder, f.name));
    % Be careful here, if top row is empty in csv, it often throws errors
    % Region filter
    if ~ismember('Region', T.Properties.VariableNames)
        warning('%s: no "Region" column → skipping', f.name);
        continue
    end
    T = T(strcmpi(T.Region, REGION), :);
    if isempty(T)
        continue
    end
    
    % Cell classification by Mike Muniak and Mao Lab also returns their own
    % boolean inside mask column (see github link). We chose not to use
    % this.
    % insideMask = T.inside_mask;
    % if ~islogical(insideMask)
    %     insideMask = logical(insideMask);   
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
    allData = [allData; allDataByType.(ctKeys{k})];
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

%% 2a. Build Voxels
% Convert to voxel‐indices and define integer bin edges
% Per-type voxel indices
voxelIdxByType = struct();
for k = 1:numel(ctKeys)
    ct = ctKeys{k};
    coords_um = allDataByType.(ct);
    if isempty(coords_um)
        voxelIdxByType.(ct) = zeros(0,3);
    else
        voxelIdxByType.(ct) = round(coords_um / voxelSize);  % [Xi Yi Zi]
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

% Per-type accumulators (for bootstrap running sums/means/vars)
accumByType = struct();
for k = 1:numel(ctKeys)
    ct = ctKeys{k};
    accumByType.(ct) = struct( ...
        'runningSum',           zeros(gridSz), ...
        'runningSum_unsmooth',  zeros(gridSz));
end

% summary
fprintf('Grid size: [%d %d %d]\n', gridSz(1), gridSz(2), gridSz(3));
for k = 1:numel(ctKeys)
    ct = ctKeys{k};
    fprintf('%s: %d points, voxelIdx size %dx3\n', ct, size(allDataByType.(ct),1), size(voxelIdxByType.(ct),1));
end

%% 2b. Exclude voxels outside caudoputamen
% Build a striatum mask in CCFv3 (25 µm) and map to your 150 µm grid
% Start with the local 150-µm grid (gridSz) and shift it by imin to get global 150-µm indices,
% then convert those indices to physical positions in microns by multiplying by voxelSize.
% Map these um coordinates to the 25-µm atlas by taking round(…/25) (nearest neighbor).
% Next, use ndgrid/sub2ind to resample the 25-µm striatum mask (`M25`) onto your 150-µm grid—this yields `mask150`, which flags whether each analysis voxel lies inside striatum.
% Optionally, apply a core-only step (L∞ erosion) so that a full 150-um cube centered at any kept voxel is guaranteed to be entirely within tissue.
%M25 = the actual 25 µm striatum mask (3-D logical). ml_25/dv_25/ap_25 = index lookups into M25. mask150 = the mask on your 150 µm grid after resampling M25.

cpFile = 'structure_672.nrrd';  % 0 = outside, 1 = inside
% Load NRRD and make a logical mask at 25 µm ---
V = nrrdread(cpFile);          % V: [z y x] = [AP DV ML] - correct with dimensions and Allen API

% Reorder to (ML, DV, AP) to match the rest of your code
M25 = permute(V > 0, [3 2 1]); % now size(M25) = [ML, DV, AP]

assert(isequal(size(V), [528 320 456]), 'Unexpected V dims; check nrrdread.');
assert(isequal(size(M25), [456 320 528]), 'Unexpected M25 dims; check permute/order.');

% grab dims if you want them explicitly
ML = size(M25,1);
DV = size(M25,2);
AP = size(M25,3);

% Build bin-center coordinates (um) for the shared 150 µm grid
% Use the global grid defined earlier (from all types combined)
nX = gridSz(1); nY = gridSz(2); nZ = gridSz(3);

ml_idx150 = imin(1) : imax(1);
dv_idx150 = imin(2) : imax(2);
ap_idx150 = imin(3) : imax(3);

ml_um = ml_idx150 * voxelSize;
dv_um = dv_idx150 * voxelSize;
ap_um = ap_idx150 * voxelSize;

% index, MATLAB starts index counting by +1. NRRD coordinates are 0-based voxel centers; MATLAB arrays are 1-based.
ml_25 = round(ml_um / 25) + 1;
dv_25 = round(dv_um / 25) + 1;
ap_25 = round(ap_um / 25) + 1;

assert(all(ml_25>=1 & ml_25<=size(M25,1)), 'ML OOB: [%g, %g]', min(ml_25), max(ml_25));
assert(all(dv_25>=1 & dv_25<=size(M25,2)), 'DV OOB: [%g, %g]', min(dv_25), max(dv_25));
assert(all(ap_25>=1 & ap_25<=size(M25,3)), 'AP OOB: [%g, %g]', min(ap_25), max(ap_25));

fprintf('ML_um range: [%.1f, %.1f] | atlas max center: %.1f\n', ...
    min(ml_um), max(ml_um), (size(M25,1)-1)*25);

% Sample the 25 µm mask onto the 150 µm grid (nearest-neighbor)
[IX, IY, IZ] = ndgrid(ml_25, dv_25, ap_25);  % order matches M25: (ML, DV, AP)
mask150 = M25(sub2ind(size(M25), IX, IY, IZ));  % [nX nY nZ] logical

% Core-only mask
rim_um   = 50;                 % 50 µm "inset" from boundary
rim_vox  = ceil(rim_um / 25);  % = 2 voxels at 25 µm

d_inf = bwdist(~M25, 'chessboard');   % L∞ distance in 25-µm voxels
core25 = d_inf >= rim_vox;            % at least rim_um from boundary

mask150_core = core25(sub2ind(size(M25), IX, IY, IZ));   % [Xn Yn Zn] type logical

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
dvent_inf = bwdist(vent25, 'chessboard');   % 0 inside c

ventRim_um  = 175;                           % desired ventricle buffer (µm)
ventRim_vox = ceil(ventRim_um / 25);        

% Keep only voxels at least 175 µm away from ventricles
ventFar25 = dvent_inf >= ventRim_vox;       % 1 = ≥ 175 µm from ventricle

% Resample ventricle buffer to 150-µm grid using same IX,IY,IZ map
ventFar150 = ventFar25(sub2ind(size(ventFar25), IX, IY, IZ));  % [Xx Yn Zn] logical

% Combine CP core mask with ventricle buffer
mask150 = mask150 & ventFar150; 


%% 3. Map cells to voxels: Bootstrap → Bin → Smooth (per cell type)
% per-type outputs
px_bootByType = struct();   % [nX x Nboot]
py_bootByType = struct();   % [nY x Nboot]
pz_bootByType = struct();   % [nZ x Nboot]

meanVol_unsmoothByType = struct();
meanVolByType          = struct();

for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    H  = hemiDataByType.(ct);
    M  = numel(H);

    % Guard: accumulators should be zeroed in earlier block; fail fast if not
    if any(accumByType.(ct).runningSum(:))      || ...
            any(accumByType.(ct).runningSum_unsmooth(:))
        error('accumByType.%s accumulators are non-zero. Clear or re-create them before bootstrapping.', ct);
    end

    % Allocate per-type bootstrap profiles
    px_bootByType.(ct) = zeros(nX, Nboot);
    py_bootByType.(ct) = zeros(nY, Nboot);
    pz_bootByType.(ct) = zeros(nZ, Nboot);

    for b = 1:Nboot
        % Resample hemisphere buckets with replacement
        if M == 0
            samplePts = [];  % no data for this type
        else
            %pick = 1:M to try the non-bootstrapped
            pick = randsample(M, M, true);
            samplePts = zeros(0,3);
            for k = 1:M
                C = H(pick(k)).coords;                 % [ML DV AP] (µm), unflipped
                if isempty(C), continue; end
                if H(pick(k)).hemi == 'L'              % flip left → right
                    C(:,1) = -C(:,1) + width;          % ML' = -ML + 2*midline
                end
                samplePts = [samplePts; C];
            end
        end

        % Bin to 3D histogram on the shared grid
        if isempty(samplePts)
            counts = zeros(gridSz);
        else
            sampIdx = round(samplePts / voxelSize);    
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
                error('Unexpected non-zeros on last-edge planes: [sx=%d sy=%d sz=%d]. Check pts / rounding / imin/imax.', sx, sy, sz);
            end
        end

        % Final safety assert: counts now matches your gridSz
        assert(isequal(size(counts), gridSz), 'counts size %s != gridSz %s', mat2str(size(counts)), mat2str(gridSz));

        % Calculate DENSITY, NOT PROBABILITY
        tot = sum(counts(:));
        if tot == 0
            pVol_unsmooth = zeros(gridSz);
        else
            pVol_unsmooth = (counts./M)./voxelVol_mm3 ;% ./tot calculation would give a probability, which we do not want here
        end

        % % 1D profiles from UNSMOOTHED density volume (mask outside as
        % NaN), commented out because we are using smoothed densities)
        % P_b = pVol_unsmooth;
        % P_b(~mask150) = NaN;
        %
        % px_bootByType.(ct)(:, b) = squeeze( mean(P_b, [2 3], 'omitnan') );  % ML
        % py_bootByType.(ct)(:, b) = squeeze( mean(P_b, [1 3], 'omitnan') );  % DV
        % pz_bootByType.(ct)(:, b) = squeeze( mean(P_b, [1 2], 'omitnan') );  % AP

        % Update UNSMOOTHED accumulators
        % I mask these values next block
        accumByType.(ct).runningSum_unsmooth = accumByType.(ct).runningSum_unsmooth + pVol_unsmooth;

        %  SMOOTHED probability volume (Gaussian sigma in voxel units)
        if tot == 0
            pVol_smooth = zeros(gridSz);
        else
            pVol_smooth = imgaussfilt3(pVol_unsmooth, filtSigma);
        end
        %  1D profiles from SMOOTHED density volume (mask outside as NaN)
        P_b = pVol_smooth;
        P_b(~mask150) = NaN;
        %
        px_bootByType.(ct)(:, b) = squeeze( mean(P_b, [2 3], 'omitnan') );  % ML
        py_bootByType.(ct)(:, b) = squeeze( mean(P_b, [1 3], 'omitnan') );  % DV
        pz_bootByType.(ct)(:, b) = squeeze( mean(P_b, [1 2], 'omitnan') );  % AP


        % Update SMOOTHED accumulators
        accumByType.(ct).runningSum = accumByType.(ct).runningSum + pVol_smooth;
    end
    % After all bootsrapping for this cell type
    % Per-type mean and std volumes (unsmoothed + smoothed)
    meanVol_unsmoothByType.(ct) = accumByType.(ct).runningSum_unsmooth / Nboot;

    meanVolByType.(ct) = accumByType.(ct).runningSum / Nboot;

end

%% 3b) OPTIONAL: Sex-split bootstrap (non-destructive; writes new variables)
doSexSplit = false;   % flip false to skip

if doSexSplit
    sexKeys = {'M','F'};

    px_bootByTypeSex = struct();
    py_bootByTypeSex = struct();
    pz_bootByTypeSex = struct();

    meanVol_unsmoothByTypeSex = struct();
    meanVolByTypeSex          = struct();

    for c = 1:numel(ctKeys)
        ct = ctKeys{c};
        H  = hemiDataByType.(ct);

        % init
        px_bootByTypeSex.(ct) = struct();  py_bootByTypeSex.(ct) = struct();  pz_bootByTypeSex.(ct) = struct();
        meanVol_unsmoothByTypeSex.(ct) = struct();  meanVolByTypeSex.(ct) = struct();

        for s = 1:numel(sexKeys)
            sx = sexKeys{s};

            % Filter hemi buckets by sex (we already have sex per bucket)
            if isempty(H)
                Hsx = H;
            else
                Hsx = H(strcmp({H.sex}, sx));
            end
            Msx = numel(Hsx);

            % Fresh accumulators per ct x sex
            runningSum          = zeros(gridSz);
            runningSum_unsmooth = zeros(gridSz);

            % Boot profiles
            px_bootByTypeSex.(ct).(sx) = zeros(nX, Nboot);
            py_bootByTypeSex.(ct).(sx) = zeros(nY, Nboot);
            pz_bootByTypeSex.(ct).(sx) = zeros(nZ, Nboot);

            for b = 1:Nboot
                % Resample sex-specific buckets
                if Msx == 0
                    samplePts = [];
                else
                    pick = randsample(Msx, Msx, true);
                    samplePts = zeros(0,3);

                    for k = 1:Msx
                        C = Hsx(pick(k)).coords;  % [ML DV AP] um, unflipped
                        if isempty(C), continue; end
                        if Hsx(pick(k)).hemi == 'L'
                            C(:,1) = -C(:,1) + width;   % flip L->R
                        end
                        samplePts = [samplePts; C];
                    end
                end

                % Bin on shared grid
                if isempty(samplePts)
                    counts = zeros(gridSz);
                else
                    sampIdx = round(samplePts / voxelSize);
                    counts  = histcnd(sampIdx, pts);
                end

                % Harmonize histcnd
                edgesSz = [numel(pts{1}), numel(pts{2}), numel(pts{3})];
                if isequal(size(counts), edgesSz)
                    sx_last = sum(counts(end,   :,   :), 'all');
                    sy_last = sum(counts(:,   end,   :), 'all');
                    sz_last = sum(counts(:,     :, end), 'all');
                    if sx_last==0 && sy_last==0 && sz_last==0
                        counts = counts(1:end-1, 1:end-1, 1:end-1);
                    else
                        error('Non-zeros on last-edge planes for %s %s (b=%d).', ct, sx, b);
                    end
                end
                assert(isequal(size(counts), gridSz), 'counts size mismatch');

                % Density: normalize by sex-specific bucket count Msx (per-bucket density)
                if Msx == 0
                    pVol_unsmooth = zeros(gridSz);
                else
                    pVol_unsmooth = (counts ./ Msx) ./ voxelVol_mm3;
                end

                runningSum_unsmooth = runningSum_unsmooth + pVol_unsmooth;

                % Smooth
                if Msx == 0
                    pVol_smooth = zeros(gridSz);
                else
                    pVol_smooth = imgaussfilt3(pVol_unsmooth, filtSigma);
                end

                % 1D profiles from smoothed + masked
                P_b = pVol_smooth;
                P_b(~mask150) = NaN;

                px_bootByTypeSex.(ct).(sx)(:, b) = squeeze(mean(P_b, [2 3], 'omitnan'));
                py_bootByTypeSex.(ct).(sx)(:, b) = squeeze(mean(P_b, [1 3], 'omitnan'));
                pz_bootByTypeSex.(ct).(sx)(:, b) = squeeze(mean(P_b, [1 2], 'omitnan'));

                runningSum = runningSum + pVol_smooth;
            end

            meanVol_unsmoothByTypeSex.(ct).(sx) = runningSum_unsmooth / Nboot;
            meanVolByTypeSex.(ct).(sx)          = runningSum / Nboot;

            fprintf('%s %s: %d hemi buckets\n', ct, sx, Msx);
        end
    end
end


%% 4) Expected number of cells per hemisphere within grid

for c = 1:numel(ctKeys)
    ct = ctKeys{c};

    % summing meanVol gives the
    s_uns_pre = sum(meanVol_unsmoothByType.(ct)(:)) * voxelVol_mm3;
    s_smo_pre = sum(meanVolByType.(ct)(:)) * voxelVol_mm3;
    fprintf('[%s] Pre-mask:  sum(meanVol_unsmooth)=%.6f,  sum(meanVol)=%.6f\n', ...
        ct, s_uns_pre, s_smo_pre);

    % Apply the striatum mask (use NaN so masked voxels don’t affect means/plots with omitnan)
    mv_uns = meanVol_unsmoothByType.(ct);
    mv_smo = meanVolByType.(ct);

    mv_uns(~mask150) = NaN;
    mv_smo(~mask150) = NaN;

    % Write back
    meanVol_unsmoothByType.(ct) = mv_uns;
    meanVolByType.(ct)          = mv_smo;

    % AFTER masking with NaNs (conditional mass inside mask)
    s_uns_post = nansum(mv_uns(:)) * voxelVol_mm3;
    s_smo_post = nansum(mv_smo(:)) * voxelVol_mm3;
    fprintf('[%s] Post-mask: nansum(meanVol_unsmooth)=%.6f, nansum(meanVol)=%.6f\n', ...
        ct, s_uns_post, s_smo_post);
end

% If you want to mask raw integer counts per bootstrap too, do it inside the bootstrap loop
% right after `counts = ...` (grids are aligned):
% counts(~mask150) = 0;

%% 5) Sanity checks & quick visuals (per cell type)

% Pick AP slice to inspect: use BregmaIdx if available; else nearest to Bregma (µm)
if exist('BregmaIdx','var')
    sliceIdx = BregmaIdx;
else
    % ap_idx150 are the global 150-µm voxel indices along AP (built earlier in mask block)
    targetAPidx = round(Bregma / voxelSize);
    [~, sliceIdx] = min(abs(ap_idx150 - targetAPidx));
end

% 0) One-time assertion: mask must match grid & any one meanVol
refCT = ctKeys{find(~cellfun(@isempty, struct2cell(meanVolByType)), 1, 'first')};
assert(~isempty(refCT), 'No meanVol found in meanVolByType.');
assert(isequal(size(mask150), size(meanVolByType.(refCT))), 'mask150 and meanVol size mismatch');

% 1) View the mask alone (same for all CTs)
figure('Name','Striatum mask (mask150) — coronal slice');
imagesc(squeeze(mask150(:,:,sliceIdx))');   % LM × DV view
set(gca,'XDir','reverse'); axis image off;
colormap(gray);
title(sprintf('mask150 at AP slice %d', sliceIdx));

% 2) Overlay mask outline on each CT’s mean volume for alignment check
for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    MV = meanVolByType.(ct);
    if isempty(MV), continue; end  % skip empty types

    figure('Name', sprintf('[%s] MeanVol with striatum outline', ct));
    img = squeeze(MV(:,:,sliceIdx))';
    imagesc(img, 'AlphaData', ~isnan(img));    % hide NaNs outside mask if present
    set(gca,'XDir','reverse'); axis image off; colormap(flipud(gray(256)));
    hold on;
    edge = bwperim(squeeze(mask150(:,:,sliceIdx))');
    [hY,hX] = find(edge);
    plot(hX, hY, '.', 'MarkerSize', 4);
    title(sprintf('[%s] AP slice %d with striatum outline', ct, sliceIdx));
end

% 3) Quick coverage stats for this slice (mask-only; same for all CTs)
insideFrac = nnz(squeeze(mask150(:,:,sliceIdx))) / numel(squeeze(mask150(:,:,sliceIdx)));
fprintf('Slice %d: %0.1f%% of bins are inside the striatum mask.\n', sliceIdx, insideFrac*100);


%% 6) 1D density profiles across axes (per cell type) with shaded SE

% Axes in mm (shared for all CTs; built from previously computed bin centers)
x_mm = (ml_um - midline) / 1000;      % ML: mm from midline
y_mm = (dv_um - dv_um(1)) / 1000;     % DV: mm from dorsal-most plane
z_mm = -(ap_um - Bregma) / 1000;      % AP: mm from Bregma (anterior positive)
% These indecies are aligned with the densities because ml_um is defined
% earlier using ml_idx150, which is using nX (from gridSz)
% Outputs (stored for later use, e.g., 6b)
densityXByType = struct();  semXByType = struct();
densityYByType = struct();  semYByType = struct();
densityZByType = struct();  semZByType = struct();

for c = 1:numel(ctKeys)
    ct = ctKeys{c};

    pxB = px_bootByType.(ct);  % [nX x Nboot]
    pyB = py_bootByType.(ct);  % [nY x Nboot]
    pzB = pz_bootByType.(ct);  % [nZ x Nboot]

    % Skip if no bootstrap profiles for this type
    if isempty(pxB) && isempty(pyB) && isempty(pzB)
        continue;
    end

    % Bootstrap means (per bin)
    densityXByType.(ct) = nanmean(pxB, 2);
    densityYByType.(ct) = nanmean(pyB, 2);
    densityZByType.(ct) = nanmean(pzB, 2);

    % Bootstrap Standard errors of bootstrapped mean estimates (SD across bootstrap reps)
    semXByType.(ct) = nanstd(pxB, [], 2);
    semYByType.(ct) = nanstd(pyB, [], 2);
    semZByType.(ct) = nanstd(pzB, [], 2);

    % Common colors for this CT's profiles
    baseColor   = [0 0.4470 0.7410];                  % mean line + underfill
    lightFactor = 0.5;                                % blend toward white
    semColor    = baseColor + (1 - baseColor)*lightFactor;

    % ML
    x_full  = x_mm(:);
    y_full  = densityXByType.(ct)(:);
    se_full = semXByType.(ct)(:);

    assert(numel(x_full) == numel(y_full) && numel(y_full) == numel(se_full), ...
        '[%s ML] Length mismatch between x, mean, and SE.', ct);

    % Drop NaNs/Infs for plotting only (outside-mask bins)
    valid = isfinite(y_full) & isfinite(se_full);
    x = x_full(valid);
    y = y_full(valid);
    se = se_full(valid);

    if ~isempty(x)
        upper = y + se;   % mean + SE


        figure('Name', sprintf('[%s] Density — ML', ct)); hold on; box on; grid on;

        % 1) Fill under mean (0 → mean)
        xx = [x; flipud(x)]; % mm coordinates, the concatenation allows for the filling
        yy = [y; zeros(size(y))]; % densities
        fill(xx, yy, baseColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

        % 2) Fill standard error band (mean → mean+SE)
        xx2 = [x; flipud(x)];
        yy2 = [upper; flipud(y)];
        fill(xx2, yy2, semColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        % 3) Lines
        plot(x, y,     'Color', baseColor, 'LineWidth', 2);    % mean
        plot(x, upper, 'Color', semColor,  'LineWidth', 1.2);  % mean + SE

        ylim([0, 1750]);

        xlabel('ML (mm from midline)');
        ylabel('Detection density');
        title(sprintf('[%s] ML density with SE band', ct));
    end

    % DV
    x_full  = y_mm(:);
    y_full  = densityYByType.(ct)(:);
    se_full = semYByType.(ct)(:);

    assert(numel(x_full) == numel(y_full) && numel(y_full) == numel(se_full), ...
        '[%s DV] Length mismatch between x, mean, and SE.', ct);

    valid = isfinite(y_full) & isfinite(se_full);
    x = x_full(valid);
    y = y_full(valid);
    se = se_full(valid);

    if ~isempty(x)
        upper = y + se;

        figure('Name', sprintf('[%s] Density — DV', ct)); hold on; box on; grid on;

        xx = [x; flipud(x)];
        yy = [y; zeros(size(y))];
        fill(xx, yy, baseColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

        xx2 = [x; flipud(x)];
        yy2 = [upper; flipud(y)];
        fill(xx2, yy2, semColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        plot(x, y,     'Color', baseColor, 'LineWidth', 2);
        plot(x, upper, 'Color', semColor,  'LineWidth', 1.2);

        ylim([0, 2500]);
        xlabel('DV (mm from dorsal)');
        ylabel('Detection density');
        title(sprintf('[%s] DV density with SE band', ct));
    end

    % AP
    x_full  = z_mm(:);
    y_full  = densityZByType.(ct)(:);
    se_full = semZByType.(ct)(:);

    assert(numel(x_full) == numel(y_full) && numel(y_full) == numel(se_full), ...
        '[%s AP] Length mismatch between x, mean, and SE.', ct);

    valid = isfinite(y_full) & isfinite(se_full);
    x = x_full(valid);
    y = y_full(valid);
    se = se_full(valid);

    if ~isempty(x)
        upper = y + se;


        figure('Name', sprintf('[%s] Density — AP', ct)); hold on; box on; grid on;

        xx = [x; flipud(x)];
        yy = [y; zeros(size(y))];
        fill(xx, yy, baseColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

        xx2 = [x; flipud(x)];
        yy2 = [upper; flipud(y)];
        fill(xx2, yy2, semColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        plot(x, y,     'Color', baseColor, 'LineWidth', 2);
        plot(x, upper, 'Color', semColor,  'LineWidth', 1.2);

        ylim([0, 2250]);
        set(gca, 'XDir','reverse');  % preserve AP convention
        xlabel('AP (mm from Bregma)');
        ylabel('Detection density');
        title(sprintf('[%s] AP density with SE band', ct));
    end
end

%% 6b Histogram of 1D density per anatomic axis with shaded error
%{
% Precompute bin widths (assumes roughly uniform spacing)
dx_mm = median(diff(x_mm));
dy_mm = median(diff(y_mm));
dz_mm = median(diff(z_mm));

for c = 1:numel(ctKeys)
    ct = ctKeys{c};

    % Require densities/SEs to exist for this CT
    if ~isfield(densityXByType, ct) || isempty(densityXByType.(ct))
        continue;
    end

    % ML axis data 
    Xx  = x_mm(:);
    Yx  = densityXByType.(ct)(:);
    SEx = semXByType.(ct)(:);

    assert(numel(Xx) == numel(Yx) && numel(Yx) == numel(SEx), ...
        '[%s ML] Length mismatch between axis, mean, SE.', ct);


    Yx(~isfinite(Yx))   = 0;
    SEx(~isfinite(SEx)) = 0;

    % DV axis data 
    Xy  = y_mm(:);
    Yy  = densityYByType.(ct)(:);
    SEy = semYByType.(ct)(:);
    assert(numel(Xy) == numel(Yy) && numel(Yy) == numel(SEy), ...
        '[%s ML] Length mismatch between axis, mean, SE.', ct);



    Yy(~isfinite(Yy))   = 0;
    SEy(~isfinite(SEy)) = 0;

    % AP axis data 
    Xz  = z_mm(:);
    Yz  = densityZByType.(ct)(:);
    SEz = semZByType.(ct)(:);

    assert(numel(Xz) == numel(Yz) && numel(Yz) == numel(SEz), ...
        '[%s ML] Length mismatch between axis, mean, SE.', ct);

    Yz(~isfinite(Yz))   = 0;
    SEz(~isfinite(SEz)) = 0;


    % Figure layout 
    fig = figure('Name', sprintf('[%s] Marginal densities (histogram style)', ct));
    t = tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');
    sgtitle(t, sprintf('%s — Cell Densities', ct), ...
        'FontSize', 12, 'FontWeight','bold');

    % ML
    ax1 = nexttile; hold(ax1,'on');

    % Bars = mean density
    bh1 = bar(ax1, Xx, Yx, 1.0);
    setBarProps(bh1);

    % Base bar color (handle default / flat cases)
    if numel(bh1) == 1 && ~isempty(bh1.FaceColor) && ~isequal(bh1.FaceColor, "flat")
        baseColor = bh1.FaceColor;
    else
        baseColor = [0 0.4470 0.7410];  % MATLAB default blue
    end

    % Lighter shade for SE cap
    alphaShade = 0.5;  % 0=original, 1=white
    shadeColor = baseColor + (1 - baseColor) * alphaShade;

    % Width of shaded cap inside each bin
    w = dx_mm * 1;

    % SE caps: [Y, Y+SE] above each bar
    for i = 1:numel(Xx)
        if SEx(i) > 0 && Yx(i) > 0
            xL = Xx(i) - w/2;
            xR = Xx(i) + w/2;
            yB = Yx(i);
            yT = Yx(i) + SEx(i);

            patch(ax1, [xL xR xR xL], [yB yB yT yT], shadeColor, ...
                'EdgeColor','none', 'FaceAlpha',0.4);
        end
    end

    xlim(ax1, [Xx(1)-dx_mm/2, Xx(end)+dx_mm/2]);
    grid(ax1,'on'); box(ax1,'off');
    xlabel(ax1,'ML (mm from midline)');
    ylabel(ax1,'Detection density');
    title(ax1,'ML');

    %DV
    ax2 = nexttile; hold(ax2,'on');

    % Bars = mean density
    bh2 = bar(ax2, Xy, Yy, 1.0);
    setBarProps(bh2);

    % Match shade color to this axis's bar color (in case color order changes)
    if numel(bh2) == 1 && ~isempty(bh2.FaceColor) && ~isequal(bh2.FaceColor, "flat")
        baseColor2 = bh2.FaceColor;
    else
        baseColor2 = baseColor;
    end
    shadeColor2 = baseColor2 + (1 - baseColor2) * alphaShade;

    w2 = dy_mm * 1;

    for i = 1:numel(Xy)
        if SEy(i) > 0 && Yy(i) > 0
            xL = Xy(i) - w2/2;
            xR = Xy(i) + w2/2;
            yB = Yy(i);
            yT = Yy(i) + SEy(i);

            patch(ax2, [xL xR xR xL], [yB yB yT yT], shadeColor2, ...
                'EdgeColor','none', 'FaceAlpha',0.4);
        end
    end

    xlim(ax2, [Xy(1)-dy_mm/2, Xy(end)+dy_mm/2]);
    grid(ax2,'on'); box(ax2,'off');
    xlabel(ax2,'DV (mm from dorsal)');
    ylabel(ax2,'Detection density');
    title(ax2,'DV');

    %AP
    ax3 = nexttile; hold(ax3,'on');

    % Bars = mean density
    bh3 = bar(ax3, Xz, Yz, 1.0);
    setBarProps(bh3);
    set(ax3,'XDir','reverse');

    if numel(bh3) == 1 && ~isempty(bh3.FaceColor) && ~isequal(bh3.FaceColor, "flat")
        baseColor3 = bh3.FaceColor;
    else
        baseColor3 = baseColor;
    end
    shadeColor3 = baseColor3 + (1 - baseColor3) * alphaShade;

    w3 = dz_mm * 1;

    for i = 1:numel(Xz)
        if SEz(i) > 0 && Yz(i) > 0
            xL = Xz(i) - w3/2;
            xR = Xz(i) + w3/2;
            yB = Yz(i);
            yT = Yz(i) + SEz(i);

            patch(ax3, [xL xR xR xL], [yB yB yT yT], shadeColor3, ...
                'EdgeColor','none', 'FaceAlpha',0.4);
        end
    end

    % AP axis limits (respect reversed XDir)
    Xz_f = Xz(isfinite(Xz));
    if isempty(Xz_f), Xz_f = Xz; end
    xlo  = min(Xz_f) - abs(dz_mm)/2;
    xhi  = max(Xz_f) + abs(dz_mm)/2;
    xlim(ax3, [xlo, xhi]);

    grid(ax3,'on'); box(ax3,'off');
    xlabel(ax3,'AP (mm from Bregma)');
    ylabel(ax3,'Detection density');
    title(ax3,'AP');
end
%}

%% 7) [not used] 2D cell density plots: every other AP section (2,4,6,...)
%{
for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    MV = meanVolByType.(ct);
    if isempty(MV), continue; end

    nSlices = size(MV, 3);
    evenIdx = 2:2:nSlices;
    nEven   = numel(evenIdx);
    if nEven == 0, continue; end

    % roughly square grid
    nCols = ceil(sqrt(nEven));
    nRows = ceil(nEven / nCols);

    % choose Bregma slice index
    if exist('BregmaIdx','var') && ~isempty(BregmaIdx)
        bregmaSlice = BregmaIdx;
    else
        % uses ap_idx150 from earlier mask block
        targetAPidx = round(Bregma / voxelSize);
        [~, bregmaSlice] = min(abs(ap_idx150 - targetAPidx));
    end

    % consistent color limits across panels
    if exist('clim','var') && ~isempty(clim)
        cax = clim;
    else
        mx  = max(MV, [], 'all', 'omitnan');
        cax = [0, mx];
    end

    fig = figure('Name', sprintf('%s Distribution — every other AP slice', ct));
    for k = 1:nEven
        sliceIdx = evenIdx(k);
        subplot(nRows, nCols, k);
        img = squeeze(MV(:,:,sliceIdx))';
        imagesc(img, cax);                 % shared clim
        set(gca,'XDir','reverse');         % display-only L-R flip
        axis image off;

        % mm relative to Bregma for this slice (2 * sliceSpacing because we step by 2)
        z_mm = -1 * (sliceIdx - bregmaSlice) * sliceSpacing;
        title(sprintf('AP %d  (Z = %+0.2f mm)', sliceIdx, z_mm), 'FontSize', 8);
    end
    colormap(fig, flipud(gray(256)));

    % shared colorbar
    hcb = colorbar('Position', [0.92 0.1 0.02 0.8]);
    hcb.Limits = cax;
    ylabel(hcb, sprintf('%s — Density of detected cell in voxel', ct));

    % note: every-other slice → 2 * voxelSize spacing
    sgtitle(fig, sprintf('%s Distribution — (%.0f µm steps)', ct, 2*voxelSize), ...
        'FontSize', 12, 'FontWeight', 'bold');
end
%}
%% 8) 2D cell density plots for all IN subtypes every-other AP section, pooled N-quantile map
% Set how many quantiles to plot
nQuant = 15;  % <-- change this to 3,5,10,... as needed

% Pool all meanVol values (inside mask, finite) to get global quantiles & cmax
vals = [];
for c = 1:numel(ctKeys)
    ct = ctKeys{c};
    V  = meanVolByType.(ct);
    if isempty(V), continue; end
    vals = [vals; V(mask150 & isfinite(V))]; %
end
assert(~isempty(vals), 'No pooled values found to compute quantiles.');

qEdges = quantile(vals, linspace(0,1,nQuant+1));   % pooled N-quantile edges
cmax   = max(vals);
clim   = [0, cmax];

% Reference size / slices
refCT   = ctKeys{find(~cellfun(@isempty, struct2cell(meanVolByType)), 1, 'first')};
nSlices = size(meanVolByType.(refCT), 3);
evenIdx = 2:2:nSlices;
nEven   = numel(evenIdx);

% Drop the very last section/slice
if nEven > 0
    evenIdx = evenIdx(1:end-1);
    nEven   = numel(evenIdx);
end
if nEven == 0
    warning('No slices to plot after dropping the last slice.');
    return
end

% Choose Bregma-based slice index (if not already defined)
if ~exist('BregmaIdx','var') || isempty(BregmaIdx)
    targetAPidx = round(Bregma / voxelSize);
    [~, BregmaIdx] = min(abs(ap_idx150 - targetAPidx));
end

% Tiled layout: ROWS = cell types, COLS = slices  <<< HORIZONTAL >>>
t = tiledlayout(numel(ctKeys), nEven, 'TileSpacing','compact', 'Padding','loose');
sgtitle(t, sprintf('%d-quantile density (300 micron steps)', nQuant), 'FontSize', 12, 'FontWeight','bold');

% Discrete N-level grayscale (lowest→highest), similar to your prior 4-bin map
g = linspace(1, 0, nQuant)';      % light→dark, can change the contrast range
cmapN = [g g g];

% Helper to bin an image into N pooled quantiles (last bin right-closed)
binToQuantN = @(x) discretize(x, qEdges);  % returns 1..nQuant (NaN outside)

useMaskFull = false;   % always use core-only mask for outline,

axFirst = [];  % for shared colorbar attachment

% Render panels (columns are slices, rows are CTs)
for r = 1:nEven                     % column = slice
    sliceIdx = evenIdx(r);
    z_mm = -1 * (sliceIdx - BregmaIdx) * sliceSpacing;  % mm from Bregma

    for c = 1:numel(ctKeys)         % row = cell type
        ct = ctKeys{c};
        MV = meanVolByType.(ct);

        tileIdx = (c-1)*nEven + r;          % linear index into layout (row-major)
        ax = nexttile(t, tileIdx);
        if isempty(axFirst), axFirst = ax; end

        if isempty(MV)
            axis(ax,'off');
            continue;
        end

        img = squeeze(MV(:,:,sliceIdx));
        img(~mask150(:,:,sliceIdx)) = NaN;  % mask outside CP

        % Bin to N pooled quantiles (NaNs stay NaN)
        cats = nan(size(img));
        valMask = isfinite(img);
        cats(valMask) = binToQuantN(img(valMask));

        % Plot categorical map (1..N), hide NaNs, reverse X for display parity
        imagesc(ax, cats', 'AlphaData', isfinite(cats'));
        set(ax,'XDir','reverse'); axis(ax,'image','off');
        colormap(ax, cmapN); caxis(ax,[0.5 nQuant+0.5]); hold(ax,'on');

        % CP outline overlay
        if useMaskFull
            BW = squeeze(mask150_precore(:,:,sliceIdx))';
        else
            BW = squeeze(mask150(:,:,sliceIdx))';
        end
        contour(ax, BW, [0.5 0.5], 'k-', 'LineWidth', 0.8);

        % Row label on leftmost column
        if r == 1
            text(ax, -0.08, 0.5, ct, ...
                'Units','normalized', ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','middle', ...
                'FontSize', 10, ...
                'FontWeight','bold', ...
                'Clipping','off');
        end

        % Column header on top row (first CT)
        if c == 1
            title(ax, sprintf('Z = %+.1f mm', z_mm), 'FontSize', 8);
        end

        hold(ax,'off');
    end
end

% Shared colorbar docked to the layout's east side
cb = colorbar(axFirst);
cb.Layout.Tile = 'east';
cb.Ticks = 1:nQuant;
cb.TickLabels = arrayfun(@(k) sprintf('Q%d',k), 1:nQuant, 'UniformOutput', false);
ylabel(cb, 'Cell Density Quantiles (pooled across all cell types)');

%% 8-sex) Atlas split by sex
%OPTIONAL: Sex-split atlas using pooled (M+F) quantiles
doPlotSexSplit = true;
% Reconstruct slice indices for sex-split plot
refCT   = ctKeys{find(~cellfun(@isempty, struct2cell(meanVolByType)), 1, 'first')};
nSlices = size(meanVolByType.(refCT), 3);

evenIdx = 2:2:nSlices;

% drop last slice to match original atlas
if numel(evenIdx) > 1
    evenIdx = evenIdx(1:end-1);
end

nEven = numel(evenIdx);
% Choose Bregma-based slice index (if not already defined)
if ~exist('BregmaIdx','var') || isempty(BregmaIdx)
    targetAPidx = round(Bregma / voxelSize);
    [~, BregmaIdx] = min(abs(ap_idx150 - targetAPidx));
end


if doPlotSexSplit && exist('meanVolByTypeSex','var')
    ctOrder  = {'SST','PV','TH'};   % desired CT order
    sexOrder = {'F','M'};          % desired sex order within each CT

    nQuant = 15;

    % pooled quantiles across BOTH sexes and all cell types
    vals = [];
    sexKeys = {'M','F'};
    for c = 1:numel(ctKeys)
        ct = ctKeys{c};
        for s = 1:numel(sexKeys)
            sx = sexKeys{s};
            if ~isfield(meanVolByTypeSex.(ct), sx), continue; end
            V = meanVolByTypeSex.(ct).(sx);
            if isempty(V), continue; end
            vals = [vals; V(mask150 & isfinite(V))];
        end
    end
    assert(~isempty(vals), 'No pooled values found for sex-split quantiles.');
    qEdges = quantile(vals, linspace(0,1,nQuant+1));

    % --- reuse your slice selection logic (refCT, evenIdx, BregmaIdx, etc.) ---
    % We assume you already computed: refCT, nSlices, evenIdx, nEven, BregmaIdx, sliceSpacing, mask150, mask150_precore.

    g = linspace(1, 0, nQuant)';  cmapN = [g g g];
    binToQuantN = @(x) discretize(x, qEdges);

    useMaskFull = false;

    % Layout: rows = sex (M/F) x cell type, cols = slices
    nRow = numel(ctOrder) * numel(sexOrder);
    t2 = tiledlayout(nRow, nEven, 'TileSpacing','compact', 'Padding','loose');

    sgtitle(t2, sprintf('Sex-split %d-quantile density (pooled F+M quantiles)', nQuant), ...
        'FontSize', 12, 'FontWeight','bold');

    axFirst = [];

    for r = 1:nEven
        sliceIdx = evenIdx(r);
        z_mm = -1 * (sliceIdx - BregmaIdx) * sliceSpacing;

        for c = 1:numel(ctOrder)
            ct = ctOrder{c};

            for s = 1:numel(sexOrder)
                sx = sexOrder{s};

                % skip gracefully if missing
                if ~isfield(meanVolByTypeSex, ct) || ~isfield(meanVolByTypeSex.(ct), sx)
                    ax = nexttile(t2, ( ( (c-1)*numel(sexOrder) + s ) - 1 )*nEven + r);
                    axis(ax,'off');
                    continue;
                end

                MV = meanVolByTypeSex.(ct).(sx);

                rowIdx  = (c-1)*numel(sexOrder) + s;   % 1..nRow, ct-major, F then M
                tileIdx = (rowIdx-1)*nEven + r;

                ax = nexttile(t2, tileIdx);
                if isempty(axFirst), axFirst = ax; end

                if isempty(MV)
                    axis(ax,'off');
                    continue;
                end

                img = squeeze(MV(:,:,sliceIdx));
                img(~mask150(:,:,sliceIdx)) = NaN;

                cats = nan(size(img));
                valMask = isfinite(img);
                cats(valMask) = binToQuantN(img(valMask));

                imagesc(ax, cats', 'AlphaData', isfinite(cats'));
                set(ax,'XDir','reverse'); axis(ax,'image','off');
                colormap(ax, cmapN); caxis(ax,[0.5 nQuant+0.5]); hold(ax,'on');

                if useMaskFull
                    BW = squeeze(mask150_precore(:,:,sliceIdx))';
                else
                    BW = squeeze(mask150(:,:,sliceIdx))';
                end
                contour(ax, BW, [0.5 0.5], 'k-', 'LineWidth', 0.8);

                % leftmost row label
                if r == 1
                    text(ax, -0.08, 0.5, sprintf('%s %s', sx, ct), ...
                        'Units','normalized', ...
                        'HorizontalAlignment','right', ...
                        'VerticalAlignment','middle', ...
                        'FontSize', 10, ...
                        'FontWeight','bold', ...
                        'Clipping','off');
                end

                % top-row titles
                if rowIdx == 1
                    title(ax, sprintf('Z = %+.1f mm', z_mm), 'FontSize', 8);
                end

                hold(ax,'off');
            end
        end
    end


cb = colorbar(axFirst);
cb.Layout.Tile = 'east';
cb.Ticks = 1:nQuant;
cb.TickLabels = arrayfun(@(k) sprintf('Q%d',k), 1:nQuant, 'UniformOutput', false);
ylabel(cb, 'Cell Density Quantiles (pooled across F+M and cell types)');
end


%% 8b) 2D Composite winner-take-all map (TH vs SST vs PV), 300 µm AP steps
% Optionally sex-disaggregated (top row = F, bottom row = M)

doPlotWTASexSplit = true;   % <-- NEW: set true for sex-split WTA, false for pooled WTA
sexOrderWTA = {'F','M'};    % top row F, bottom row M

% Cell types to include in the composite
ctPlot = ctKeys;
nCT    = numel(ctPlot);

if nCT == 0
    warning('No TH/SST/PV in ctKeys; skipping winner-take-all composite.');
else

    % ---- Determine slice count robustly depending on mode ----
    refCT = ctPlot{find(~cellfun(@isempty, struct2cell(meanVolByType)), 1, 'first')};
    if isempty(refCT)
        % If pooled meanVolByType is empty but sex struct exists, try to infer from it
        if exist('meanVolByTypeSex','var') && isfield(meanVolByTypeSex, ctPlot{1}) ...
                && isfield(meanVolByTypeSex.(ctPlot{1}), sexOrderWTA{1}) ...
                && ~isempty(meanVolByTypeSex.(ctPlot{1}).(sexOrderWTA{1}))
            refCT = ctPlot{1};
            nSlices = size(meanVolByTypeSex.(refCT).(sexOrderWTA{1}), 3);
        else
            warning('meanVolByType is empty and cannot infer nSlices from meanVolByTypeSex; skipping WTA composite.');
            return;
        end
    else
        nSlices = size(meanVolByType.(refCT), 3);
    end

    % Every other slice (300 µm steps) as in Block 8
    evenIdx = 2:2:nSlices;

    % Optionally drop the last slice for cleaner layout (match Block 8 style)
    if ~isempty(evenIdx)
        evenIdx = evenIdx(1:end-1);
    end
    nEven = numel(evenIdx);

    if nEven == 0
        warning('No slices selected for winner-take-all composite.');
        return;
    end

    % Ensure BregmaIdx is defined
    if ~exist('BregmaIdx','var') || isempty(BregmaIdx)
        targetAPidx = round(Bregma / voxelSize);
        [~, BregmaIdx] = min(abs(ap_idx150 - targetAPidx));
    end

    % Define distinct RGB colors for each CT
    ctColor = struct();
    ctColor.SST = [0.1333 0.2588 0.5647];
    ctColor.PV  = [0.0000 0.6157 0.5529];
    ctColor.TH  = [0.6039 0.2275 0.5882];

    % ---- Figure layout ----
    if doPlotWTASexSplit
        if ~exist('meanVolByTypeSex','var')
            warning('doPlotWTASexSplit=true but meanVolByTypeSex not found; falling back to pooled WTA.');
            doPlotWTASexSplit = false;
        end
    end

    if doPlotWTASexSplit
        figure('Name','Max Density Cell Type per voxel (Sex-split WTA)');
        t = tiledlayout(2, nEven, 'TileSpacing','compact', 'Padding','compact');
        sgtitle(t, 'Max-density cell type per voxel (TH/SST/PV), sex-disaggregated', ...
            'FontSize', 12, 'FontWeight','bold');
    else
        figure('Name','Max Density Cell Type per voxel (Pooled WTA)');
        t = tiledlayout(1, nEven, 'TileSpacing','compact', 'Padding','compact');
        sgtitle(t, 'Max-density cell type per voxel (TH/SST/PV)', ...
            'FontSize', 12, 'FontWeight','bold');
    end

    % ---- Main loop over slices (and sex if enabled) ----
    for idx = 1:nEven
        sliceIdx = evenIdx(idx);
        z_mm = -1 * (sliceIdx - BregmaIdx) * sliceSpacing;

        % Mask for this slice (DV x ML)
        BW = squeeze(mask150(:,:,sliceIdx))';
        mask_slice = BW;

        % Determine display size from first available volume
        % (use pooled or sex volume depending on mode)
        if doPlotWTASexSplit
            % Find a CT/sex with data
            firstMV = [];
            for c = 1:nCT
                ct = ctPlot{c};
                for s = 1:numel(sexOrderWTA)
                    sx = sexOrderWTA{s};
                    if isfield(meanVolByTypeSex, ct) && isfield(meanVolByTypeSex.(ct), sx) ...
                            && ~isempty(meanVolByTypeSex.(ct).(sx))
                        firstMV = meanVolByTypeSex.(ct).(sx);
                        break;
                    end
                end
                if ~isempty(firstMV), break; end
            end
            if isempty(firstMV)
                warning('No non-empty sex volumes found for slice %d; skipping.', sliceIdx);
                continue;
            end
        else
            firstMV = meanVolByType.(ctPlot{1});
            if isempty(firstMV)
                warning('Pooled meanVolByType empty for WTA; skipping.');
                break;
            end
        end

        sliceSz = [size(firstMV,2), size(firstMV,1)]; % [nY nX] after transpose (DV x ML)

        if doPlotWTASexSplit
            for s = 1:numel(sexOrderWTA)
                sx = sexOrderWTA{s};      % 'F' then 'M'
                row = s;                  % row 1 = F, row 2 = M
                tileIdx = (row-1)*nEven + idx;

                % Build M = [DV x ML x nCT] for this sex and slice
                M = nan([sliceSz, nCT]);

                for c = 1:nCT
                    ct = ctPlot{c};

                    if ~isfield(meanVolByTypeSex, ct) || ~isfield(meanVolByTypeSex.(ct), sx)
                        continue
                    end
                    Vct = meanVolByTypeSex.(ct).(sx);
                    if isempty(Vct)
                        continue
                    end

                    % Vct is [ML x DV x AP]; transpose for display to [DV x ML]
                    M(:,:,c) = squeeze(Vct(:,:,sliceIdx))';
                end

                % Winner-take-all on this sex
                mask3D = repmat(mask_slice, [1 1 nCT]);
                outside = ~mask3D;
                M(outside) = NaN;

                M2 = M;
                M2(~isfinite(M2)) = -Inf;
                [maxVal, maxIdx] = max(M2, [], 3);

                noWinner = (~mask_slice) | (~any(isfinite(M),3)) | (maxVal <= 0);

                alphaImg = ones(size(mask_slice));
                alphaImg(noWinner) = 0;

                [nY, nX, ~] = size(M);
                imgR = ones(nY, nX); imgG = ones(nY, nX); imgB = ones(nY, nX);

                for c = 1:nCT
                    ct = ctPlot{c};
                    col = ctColor.(ct);
                    mask_ct = (maxIdx == c) & (maxVal > 0);
                    imgR(mask_ct) = col(1);
                    imgG(mask_ct) = col(2);
                    imgB(mask_ct) = col(3);
                end
                imgRGB = cat(3, imgR, imgG, imgB);

                ax = nexttile(t, tileIdx);
                imagesc(ax, imgRGB, 'AlphaData', alphaImg);
                set(ax, 'XDir','reverse');
                axis(ax, 'image', 'off');
                set(ax, 'Color', 'none');

                % Titles only on top row
                if row == 1
                    title(ax, sprintf('Z = %+0.1f mm', z_mm), 'FontSize', 8);
                end

                % Row labels at first column
                if idx == 1
                    text(ax, -0.08, 0.5, sx, ...
                        'Units','normalized', ...
                        'HorizontalAlignment','right', ...
                        'VerticalAlignment','middle', ...
                        'FontSize', 10, ...
                        'FontWeight','bold', ...
                        'Clipping','off');
                end
            end

        else
            % ---- Original pooled behavior (one row) ----
            M = nan([sliceSz, nCT]);
            for c = 1:nCT
                ct  = ctPlot{c};
                Vct = meanVolByType.(ct);
                if isempty(Vct), continue; end
                M(:,:,c) = squeeze(Vct(:,:,sliceIdx))';
            end

            mask3D = repmat(mask_slice, [1 1 nCT]);
            outside = ~mask3D;
            M(outside) = NaN;

            M2 = M;
            M2(~isfinite(M2)) = -Inf;
            [maxVal, maxIdx] = max(M2, [], 3);

            noWinner = (~mask_slice) | (~any(isfinite(M),3)) | (maxVal <= 0);
            alphaImg = ones(size(mask_slice));
            alphaImg(noWinner) = 0;

            [nY, nX, ~] = size(M);
            imgR = ones(nY, nX); imgG = ones(nY, nX); imgB = ones(nY, nX);

            for c = 1:nCT
                ct = ctPlot{c};
                col = ctColor.(ct);
                mask_ct = (maxIdx == c) & (maxVal > 0);
                imgR(mask_ct) = col(1);
                imgG(mask_ct) = col(2);
                imgB(mask_ct) = col(3);
            end
            imgRGB = cat(3, imgR, imgG, imgB);

            ax = nexttile(t, idx);
            imagesc(ax, imgRGB, 'AlphaData', alphaImg);
            set(ax, 'XDir','reverse');
            axis(ax, 'image', 'off');
            set(ax, 'Color', 'none');
            title(ax, sprintf('Z = %+0.1f mm', z_mm), 'FontSize', 8);
        end
    end

    % ---- Legend-like annotation (works for both modes) ----
    axLeg = axes('Position',[0.88 0.80 0.10 0.15]);
    axis(axLeg,'off'); hold(axLeg,'on');
    y0 = 0.8; dy = 0.25;
    fields = ctPlot;
    for i = 1:numel(fields)
        ct = fields{i};
        col = ctColor.(ct);
        rectangle(axLeg, 'Position',[0.05, y0 - i*dy, 0.15, 0.24], ...
            'FaceColor', col, 'EdgeColor','none');
        text(axLeg, 0.25, y0 - i*dy + 0.25, ct, ...
            'Units','normalized', 'FontSize',8, ...
            'VerticalAlignment','middle');
    end

end


%% 9) Statistical analysis corresponding to 1D density plots (Mixed-effects regression per cell type and axis by hemisphere)

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
useSexCovariate  = true;   % include sex as fixed effect
testSexInteraction = true; % include X:sex interaction (set false for main-effects only)
alphaFDR             = 0.05;   % BH-FDR q

lmeModels = struct();     % optional
% Collect sex-interaction p-values across ct × axis for BH-FDR
statsRows = {};  % each row: {ct, axis, nHemi, raw_p_int, coefName, beta_int, SE_int}
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
        byHemi(hCount).sex    = char(H(k).sex);  % 'M' or 'F'
        sx = char(H(k).sex);
        if ~ismember(sx, {'M','F'})
            warning('[%s] hemiID %s has unexpected sex "%s"; setting to missing.', ct, byHemi(hCount).hemiID, sx);
            byHemi(hCount).sex = 'U';
        end
        byHemi(hCount).px     = px;
        byHemi(hCount).py     = py;
        byHemi(hCount).pz     = pz;
    end

    if isempty(byHemi)
        fprintf('[%s] No per-hemisphere profiles after filtering; skipping LME.\n', ct);
        continue;
    end

    A  = struct2table(byHemi);
    nH = height(A);

    % Axis vectors (in mm)
    % Recreate from earlier bin-center coords to avoid accidental edits.
    x_mm_axis = (ml_um - midline) / 1000;        % ML: mm from midline
    y_mm_axis = (dv_um - dv_um(1)) / 1000;       % DV: mm from dorsal-most plane
    z_mm_axis = -(ap_um - Bregma) / 1000;        % AP: mm from Bregma (ant +)

    % Long-format tables: one row = one bin of one hemisphere

    % ML
    Yx = []; Xx = []; IDx = strings(0,1); Sx = strings(0,1);
    for i = 1:nH
        vx = getcol(A.px, i);
        vx = vx(:);
        assert(numel(vx) == numel(x_mm_axis), ...
            '[%s ML] Length mismatch between hemisphere %d profile and x_mm_axis.', ct, i);

        Yx  = [Yx; vx];
        Xx  = [Xx; x_mm_axis(:)];
        IDx = [IDx; repmat(string(A.hemiID{i}), numel(vx), 1)];
        Sx  = [Sx;  repmat(string(A.sex{i}),    numel(vx), 1)];
    end

    Tx = table(Yx, Xx, categorical(Sx), categorical(IDx), ...
        'VariableNames', {'Y','X','sex','hemiID'});

    % DV
    Yy = []; Xy = []; IDy = strings(0,1); Sy = strings(0,1);
    for i = 1:nH
        vy = getcol(A.py, i);
        vy = vy(:);
        assert(numel(vy) == numel(y_mm_axis), ...
            '[%s DV] Length mismatch between hemisphere %d profile and y_mm_axis.', ct, i);

        Yy  = [Yy; vy];
        Xy  = [Xy; y_mm_axis(:)];
        IDy = [IDy; repmat(string(A.hemiID{i}), numel(vy), 1)];
        Sy  = [Sy;  repmat(string(A.sex{i}),    numel(vy), 1)];
    end

    Ty = table(Yy, Xy, categorical(Sy), categorical(IDy), ...
        'VariableNames', {'Y','X','sex','hemiID'});


    % AP
    Yz = []; Xz = []; IDz = strings(0,1); Sz = strings(0,1);
    for i = 1:nH
        vz = getcol(A.pz, i);
        vz = vz(:);
        assert(numel(vz) == numel(z_mm_axis), ...
            '[%s AP] Length mismatch between hemisphere %d profile and z_mm_axis.', ct, i);

        Yz  = [Yz; vz];
        Xz  = [Xz; z_mm_axis(:)];
        IDz = [IDz; repmat(string(A.hemiID{i}), numel(vz), 1)];
        Sz  = [Sz;  repmat(string(A.sex{i}),    numel(vz), 1)];
    end

    Tz = table(Yz, Xz, categorical(Sz), categorical(IDz), ...
        'VariableNames', {'Y','X','sex','hemiID'});


    % Keep raw X for plotting
    Xx_raw = Xx;
    Xy_raw = Xy;
    Xz_raw = Xz;

    % Fit LME models
    if useRandomSlope
        reTerm = '(1 + X|hemiID)';
    else
        reTerm = '(1|hemiID)';
    end

    if useSexCovariate
        if testSexInteraction
            feTerm = 'Y ~ X*sex';
            modelStr = sprintf('Y ~ X*sex + %s', reTerm);
        else
            feTerm = 'Y ~ X + sex';
            modelStr = sprintf('Y ~ X + sex + %s', reTerm);
        end
    else
        feTerm = 'Y ~ X';
        modelStr = sprintf('Y ~ X + %s', reTerm);
    end

    formula = sprintf('%s + %s', feTerm, reTerm);

    % Random slope fit
    mdlX = fitlme(Tx, formula);
    mdlY = fitlme(Ty, formula);
    mdlZ = fitlme(Tz, formula);

    % --- Collect sex interaction stats for BH-FDR (X:sex term) ---
    if useSexCovariate && testSexInteraction
        % ML
        [pInt, coefName, bInt, seInt] = extract_sex_interaction(mdlX);
        statsRows(end+1,:) = {ct, 'ML', nH, pInt, coefName, bInt, seInt};

        % DV
        [pInt, coefName, bInt, seInt] = extract_sex_interaction(mdlY);
        statsRows(end+1,:) = {ct, 'DV', nH, pInt, coefName, bInt, seInt};

        % AP
        [pInt, coefName, bInt, seInt] = extract_sex_interaction(mdlZ);
        statsRows(end+1,:) = {ct, 'AP', nH, pInt, coefName, bInt, seInt};
    end


    % Report fixed-effect slope β_X
    CX = mdlX.Coefficients; rowX = strcmp(CX.Name,'X');
    CY = mdlY.Coefficients; rowY = strcmp(CY.Name,'X');
    CZ = mdlZ.Coefficients; rowZ = strcmp(CZ.Name,'X');

    fprintf('\n[%s] Mixed-effects by hemisphere: %s\n', ct, modelStr);

    print_effects = @(mdl, axisName, nH, nRows) ...
        fprintf(' %s: %s\n', axisName, summarize_fixed(mdl, nH, nRows));

    fprintf(' ML: %s\n', summarize_fixed(mdlX, nH, height(Tx)));
    fprintf(' DV: %s\n', summarize_fixed(mdlY, nH, height(Ty)));
    fprintf(' AP: %s\n', summarize_fixed(mdlZ, nH, height(Tz)));

    % Store models
    lmeModels.(ct).ML = mdlX;
    lmeModels.(ct).DV = mdlY;
    lmeModels.(ct).AP = mdlZ;

    % Plots: data + fixed-effect line + optional spaghetti

    annFmt = '\\beta = %.3g (SE=%.3g), p = %.3g';

    % ML
    figure('Name', sprintf('[%s] LME — ML axis', ct)); hold on; box on; grid on;
    plot(Xx_raw, Yx, '.', 'MarkerSize', 6, 'Color', 0.8*[1 1 1]);  % all points

    xg = x_mm_axis(:);  % one entry per voxel/bin along ML
    fe   = fixedEffects(mdlX);      % [b0; b1]
    yhat = fe(1) + fe(2) * xg;      % same units, same grid

    plot(xg, yhat, 'LineWidth', 2);

    if showHemiSlopes
        plot_spaghetti_lines(mdlX, Tx, xg, [0.6 0.6 0.6]);
    end

    xlabel('ML (mm from midline)');
    ylabel('Density (cells / mm^3 / hemisphere)');
    title(sprintf('[%s] ML: %s', ct, modelStr));

    bx = CX.Estimate(rowX); sx = CX.SE(rowX); pxv = CX.pValue(rowX);
    text(0.03, 0.95, sprintf(annFmt, bx, sx, pxv), ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold','BackgroundColor','w','Margin',2);

    % DV
    figure('Name', sprintf('[%s] LME — DV axis', ct)); hold on; box on; grid on;
    plot(Xy_raw, Yy, '.', 'MarkerSize', 6, 'Color', 0.8*[1 1 1]);

    fe   = fixedEffects(mdlY);
    yg = y_mm_axis(:);
    yhat = fe(1) + fe(2)*yg;
    plot(yg, yhat, 'LineWidth', 2);
    if showHemiSlopes
        plot_spaghetti_lines(mdlY, Ty, yg, [0.6 0.6 0.6]);
    end

    xlabel('DV (mm from dorsal)');
    ylabel('Density (cells / mm^3 / hemisphere)');
    title(sprintf('[%s] DV: %s', ct, modelStr));

    by = CY.Estimate(rowY); sy = CY.SE(rowY); pyv = CY.pValue(rowY);
    text(0.03, 0.95, sprintf(annFmt, by, sy, pyv), ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold','BackgroundColor','w','Margin',2);

    % AP
    figure('Name', sprintf('[%s] LME — AP axis', ct)); hold on; box on; grid on;
    plot(Xz_raw, Yz, '.', 'MarkerSize', 6, 'Color', 0.8*[1 1 1]);

    fe   = fixedEffects(mdlZ);
    zg = z_mm_axis(:);
    yhat = fe(1) + fe(2)*zg;
    plot(zg, yhat, 'LineWidth', 2);
    if showHemiSlopes
        plot_spaghetti_lines(mdlZ, Tz, zg, [0.6 0.6 0.6]);
    end

    set(gca,'XDir','reverse');  % preserve your AP convention
    xlabel('AP (mm from Bregma)');
    ylabel('Density (cells / mm^3 / hemisphere)');
    title(sprintf('[%s] AP: %s', ct, modelStr));

    bz = CZ.Estimate(rowZ); szv = CZ.SE(rowZ); pzv = CZ.pValue(rowZ);
    text(0.03, 0.95, sprintf(annFmt, bz, szv, pzv), ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',10,'FontWeight','bold','BackgroundColor','w','Margin',2);
end
% --- BH-FDR correction across all ct×axis sex-interaction tests ---
if useSexCovariate && testSexInteraction && ~isempty(statsRows)
    S = cell2table(statsRows, ...
        'VariableNames', {'cellType','axis','nHemi','p_int_raw','coefName','beta_int','se_int'});

    % Some models might not contain the term (shouldn't happen if formula is X*sex),
    % but keep robust: exclude NaNs from correction.
    p = S.p_int_raw;
    keep = isfinite(p);
    p_keep = p(keep);

    % BH-FDR (preferred: mafdr if available)
    % mafdr returns adjusted p-values (q-values) for BH-FDR when 'BHFDR' true
    q = nan(size(p));
    try
        q_keep = mafdr(p_keep, 'BHFDR', true);
        q(keep) = q_keep;
    catch
        % Fallback: manual BH implementation
        q(keep) = bh_fdr_qvalues(p_keep);
    end

    S.p_int_fdr = q;
    S.sig_fdr = S.p_int_fdr <= alphaFDR;

    % Sort for readability
    S = sortrows(S, {'p_int_fdr','p_int_raw'});

    fprintf('\n=== Sex interaction BH-FDR across ct×axis (q=%.3f), m=%d tests ===\n', ...
        alphaFDR, sum(keep));
    disp(S(:, {'cellType','axis','nHemi','p_int_raw','p_int_fdr','sig_fdr','coefName','beta_int','se_int'}));

    % Optional: stash for later reporting
    sexInteractionSummary = S;
end


% Helpers for Block 9

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

function s = summarize_fixed(mdl, nH, nRows)
C = mdl.Coefficients;

% Always report X slope
rowX = strcmp(C.Name,'X');
parts = {};
if any(rowX)
    parts{end+1} = sprintf('beta_X=%.4g (SE=%.4g), p=%.3g', ...
        C.Estimate(rowX), C.SE(rowX), C.pValue(rowX));
end

% Sex main effect: matlab typically encodes as 'sex_F' or similar depending on reference level
sexRows = startsWith(string(C.Name), "sex_");
if any(sexRows)
    for r = find(sexRows)'
        parts{end+1} = sprintf('%s=%.4g (SE=%.4g), p=%.3g', ...
            C.Name{r}, C.Estimate(r), C.SE(r), C.pValue(r));
    end
end

% Interaction term usually 'X:sex_F' etc
intRows = contains(string(C.Name), "X:sex_") | contains(string(C.Name), "sex_:X");
if any(intRows)
    for r = find(intRows)'
        parts{end+1} = sprintf('%s=%.4g (SE=%.4g), p=%.3g', ...
            C.Name{r}, C.Estimate(r), C.SE(r), C.pValue(r));
    end
end

parts{end+1} = sprintf('nHemi=%d, nRows=%d', nH, nRows);
s = strjoin(parts, '; ');
end

function [pInt, coefName, bInt, seInt] = extract_sex_interaction(mdl)
% Extract the X:sex interaction term from a fitlme model robustly.
% Returns NaN if not found.

C = mdl.Coefficients;
names = string(C.Name);

% Candidate patterns: "X:sex_*" or "sex_*:X"
isInt = contains(names, "X:sex_") | contains(names, "sex_:X");

if ~any(isInt)
    % Some MATLAB versions might encode differently; broaden slightly:
    isInt = contains(names, "X:sex") | contains(names, "sex:X");
end

if ~any(isInt)
    pInt = NaN; coefName = ""; bInt = NaN; seInt = NaN;
    return;
end

% If multiple interaction terms exist (unlikely unless sex has >2 levels),
% take the first.
idx = find(isInt, 1, 'first');

pInt = C.pValue(idx);
coefName = C.Name{idx};
bInt = C.Estimate(idx);
seInt = C.SE(idx);
end

function q = bh_fdr_qvalues(p)
% Manual Benjamini-Hochberg q-values (adjusted p-values).
% p: vector of raw p-values (finite)
p = p(:);
[ps, order] = sort(p, 'ascend');
m = numel(ps);

% BH adjusted p: q_i = p_i * m / i, then enforce monotonicity
q_sorted = ps .* (m ./ (1:m)');
q_sorted = min(q_sorted, 1);

% monotone nondecreasing when going from small to large p:
for i = m-1:-1:1
    q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
end

q = nan(size(p));
q(order) = q_sorted;
end


%% Helper functions


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
