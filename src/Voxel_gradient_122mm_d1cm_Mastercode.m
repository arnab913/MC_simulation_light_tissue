%% =========================================================
%  Diffuse Reflectance / Transmission / Ink Maps + Banana (B_wds/B_ds)
%  for ALL Conditions (3x5 ordered grid)
%
%  Assumes each .mat contains:
%    scatter_stat, t_model, sim_info
%    and (for banana) B_wds and B_ds
% ==========================================================
clear; clc; close all;

%% -------------------- USER SETTINGS ----------------------
dataFolder = "D:\MC SImulation data\vox_gradient_tissue123cmlayers";

xMin = -7.5; xMax = 7.5;
yMin = -7.5; yMax = 7.5;
zMin = 0;    zMax = 5;      % top z = 0, bottom z = 5
useLog = true;              % log-scale maps for surface/trans/ink maps

% Banana plot settings
plotBanana = true;
bananaUse = "B_wds";        % "B_wds" or "B_ds"
bananaRotate180 = true;     % apply flipud(fliplr(...))
bananaLevels = 30;
epsv = 1e-12;

% Ordering: rows and cols
rowUs1 = [30, 5, 2];              % rows: No dye, 1L dye, 2L dye
colUa5 = [0.05, 1, 7, 30, 150];    % cols: No ink, 1, 7, 30, 150
tol = 1e-9;                        % tolerance for float comparisons
% ----------------------------------------------------------

%% --------- BUILD ORDERED FILE LIST (3x5 grid) ----------
allFiles   = dir(fullfile(dataFolder, '*.mat'));
nAll       = numel(allFiles);
if nAll == 0
    error('No .mat files found in: %s', dataFolder);
end

us1 = zeros(nAll,1);
ua5 = zeros(nAll,1);

for k = 1:nAll
    S = load(fullfile(allFiles(k).folder, allFiles(k).name), 'sim_info');
    us1(k) = S.sim_info.us1;
    ua5(k) = S.sim_info.ua5;
end

nRows = numel(rowUs1);
nCols = numel(colUa5);

ordered = repmat(allFiles(1), nRows*nCols, 1);
idx = 0;

for i = 1:nRows
    for j = 1:nCols
        idx = idx + 1;

        mask = (abs(us1 - rowUs1(i)) < tol) & (abs(ua5 - colUa5(j)) < tol);
        fidx = find(mask, 1);

        if isempty(fidx)
            error('No file for us1=%g, ua5=%g', rowUs1(i), colUa5(j));
        end
        ordered(idx) = allFiles(fidx);
    end
end

files  = ordered;
nFiles = numel(files);

%% Containers
Diff_refl = cell(nFiles,1);
Trans_all = cell(nFiles,1);
Ink_all   = cell(nFiles,1);
LayerPath = cell(nFiles,1);

% Banana slices (XZ @ y≈0) stored for plotting with consistent scaling
BananaSlices = cell(nFiles,1);
BananaTitle  = strings(nFiles,1);

%% -------------------- MAIN LOOP --------------------------
for k = 1:nFiles

    % ----- Load file -----
    S = load(fullfile(files(k).folder, files(k).name));
if k==1
    fprintf('Has B_wds=%d, B_ds=%d\n', isfield(S,'B_wds'), isfield(S,'B_ds'));
end


    scatter_stat = S.scatter_stat;
    t_model      = S.t_model;

    % ----- Extract final photon data -----
    x_ot = scatter_stat(:,4);
    y_ot = scatter_stat(:,5);
    z_ot = scatter_stat(:,6);
    w_ot = scatter_stat(:,8);

    % total photons in this file
    Nphot = size(scatter_stat,1);

    % ----- Grid info from model -----
    Lx = t_model.G.Lx;
    Ly = t_model.G.Ly;
    Lz = t_model.G.Lz; %#ok<NASGU>
    Nx = t_model.G.nx;
    Ny = t_model.G.ny;

    xEdges = linspace(-Lx/2, Lx/2, Nx+1);
    yEdges = linspace(-Ly/2, Ly/2, Ny+1);

    % --------------------------------------------------
    % 1) Surface photons (diffuse reflectance)
    % --------------------------------------------------
    idx_surface = (z_ot <= zMin);

    x_surf = x_ot(idx_surface);
    y_surf = y_ot(idx_surface);
    w_surf = w_ot(idx_surface);

    n_reflect = sum(idx_surface);

    [~,~,~,idxR,idyR] = histcounts2(x_surf, y_surf, xEdges, yEdges);
    RR = zeros(Nx,Ny);
    for ii = 1:numel(idxR)
        if idxR(ii) > 0 && idyR(ii) > 0
            RR(idxR(ii), idyR(ii)) = RR(idxR(ii), idyR(ii)) + w_surf(ii);
        end
    end

    if useLog
        RR(RR <= 0) = NaN;
        RR = log10(RR);
    end

    Diff_refl{k} = struct('map', RR, 'nPhot', n_reflect, 'Nphot', Nphot);

    % --------------------------------------------------
    % 2) Transmitted photons (bottom + sides, excluding top)
    % --------------------------------------------------
    idx_trans_bottom = (z_ot >= zMax);
    idx_side_x_left  = (x_ot <= xMin);
    idx_side_x_right = (x_ot >= xMax);
    idx_side_y_front = (y_ot <= yMin);
    idx_side_y_back  = (y_ot >= yMax);

    idx_trans_sides = idx_side_x_left | idx_side_x_right | ...
                      idx_side_y_front | idx_side_y_back;

    idx_transmitted = (idx_trans_bottom | idx_trans_sides) & ~idx_surface;

    x_trans = x_ot(idx_transmitted);
    y_trans = y_ot(idx_transmitted);
    w_trans = w_ot(idx_transmitted);

    n_trans = sum(idx_transmitted);

    [~,~,~,idxT,idyT] = histcounts2(x_trans, y_trans, xEdges, yEdges);
    Trans_map = zeros(Nx,Ny);
    for ii = 1:numel(idxT)
        if idxT(ii) > 0 && idyT(ii) > 0
            Trans_map(idxT(ii), idyT(ii)) = Trans_map(idxT(ii), idyT(ii)) + w_trans(ii);
        end
    end

    if useLog
        Trans_map(Trans_map <= 0) = NaN;
        Trans_map = log10(Trans_map);
    end

    Trans_all{k} = struct('map', Trans_map, 'nPhot', n_trans, 'Nphot', Nphot);

    % --------------------------------------------------
    % 3) Photons trapped in ink (final point inside ink, not escaped)
    % --------------------------------------------------
    ink_mask = (x_ot >= -7.5) & (x_ot <= 7.5) & ...
               ((y_ot - 0).^2 + (z_ot - 1.0).^2 <= 0.3^2);
    ink_trapped_mask = ink_mask & ~idx_surface & ~idx_transmitted;

    x_ink = x_ot(ink_trapped_mask);
    y_ink = y_ot(ink_trapped_mask);
    w_ink = w_ot(ink_trapped_mask);

    n_ink = sum(ink_trapped_mask);

    [~,~,~,idxI,idyI] = histcounts2(x_ink, y_ink, xEdges, yEdges);
    Ink_map = zeros(Nx,Ny);
    for ii = 1:numel(idxI)
        if idxI(ii) > 0 && idyI(ii) > 0
            Ink_map(idxI(ii), idyI(ii)) = Ink_map(idxI(ii), idyI(ii)) + w_ink(ii);
        end
    end

    if useLog
        Ink_map(Ink_map <= 0) = NaN;
        Ink_map = log10(Ink_map);
    end

    Ink_all{k} = struct('map', Ink_map, 'nPhot', n_ink, 'Nphot', Nphot);

    % --------------------------------------------------
    % 4) Mean pathlength per layer (your existing approximation)
    % --------------------------------------------------
    z_out = scatter_stat(:,6);
    s_tot = scatter_stat(:,7);

    zL = [0.0, 0.3, 0.4, 0.5, inf];   % 3 dye layers + bulk
    nLayers = numel(zL)-1;
    Li = zeros(nLayers,1);
    Ni = zeros(nLayers,1);

    for p = 1:Nphot
        zf = z_out(p);
        s  = s_tot(p);
        if s <= 0, continue; end

        dz = zeros(nLayers,1);
        for j = 1:nLayers
            z1 = zL(j);
            z2 = zL(j+1);
            dz(j) = max(0, min(zf, z2) - z1);
        end
        dz_tot = sum(dz);

        if dz_tot > 0
            Li = Li + s * (dz / dz_tot);
            Ni = Ni + (dz > 0);
        end
    end

    Li_mean = Li ./ max(Ni,1);
    LayerPath{k} = struct('Li', Li_mean, 'zBounds', zL);

    % --------------------------------------------------
    % 5) Banana slice (voxelized fluence/pathlength) from saved B_wds / B_ds
    % --------------------------------------------------
    if plotBanana
        if isfield(S, 'sim_info') && isfield(S.sim_info, 'grid')
            g = S.sim_info.grid;
            nx=g.nx; ny=g.ny; nz=g.nz; dl=g.dl;
            Lxg=g.Lx; Lyg=g.Ly;

            xC = ((1:nx)-0.5)*dl - Lxg/2;
            yC = ((1:ny)-0.5)*dl - Lyg/2;
            zC = ((1:nz)-0.5)*dl;

            [~,iy0] = min(abs(yC - 0));

            if bananaUse == "B_wds" && isfield(S,'B_wds')
                B = S.B_wds;
            elseif bananaUse == "B_ds" && isfield(S,'B_ds')
                B = S.B_ds;
            elseif isfield(S,'B_wds')
                B = S.B_wds;
                bananaUse = "B_wds";
            elseif isfield(S,'B_ds')
                B = S.B_ds;
                bananaUse = "B_ds";
            else
                B = [];
            end

            if ~isempty(B)
                Bxz = squeeze(B(:,iy0,:))';        % nz x nx
                Bxz = log10(Bxz + epsv);
                if bananaRotate180
                    Bxz = flipud(fliplr(Bxz));
                end
                BananaSlices{k} = Bxz;
                BananaTitle(k) = files(k).name;
            end
        end
    end
end

%% ================== PLOTTING ============================
rows = nRows; cols = nCols;

% Rebuild xy grid from last loaded model
Lx = t_model.G.Lx; Ly = t_model.G.Ly;
Nx = t_model.G.nx; Ny = t_model.G.ny;
xEdges = linspace(-Lx/2, Lx/2, Nx+1);
yEdges = linspace(-Ly/2, Ly/2, Ny+1);

%% 1) Diffuse Reflectance maps
figure('Color','w','Position',[100 100 1800 900]);
for k = 1:nFiles
    subplot(rows, cols, k);
    RR = Diff_refl{k}.map;
    contourf(xEdges(1:end-1), yEdges(1:end-1), RR, 30, 'LineColor','none');
    axis equal tight; set(gca,'YDir','normal'); colorbar;
    clim([-5 5]);
    xlabel('x [cm]'); ylabel('y [cm]');

    params = parseFilename(files(k).name);
    nR = Diff_refl{k}.nPhot; Np = Diff_refl{k}.Nphot;

    title({
        sprintf('us=[%.3g %.3g %.3g]', params.us1, params.us2, params.us3)
        sprintf('ua5=%.3g us5=%.3g n5=%.2f', params.ua5, params.us5, params.n5)
        sprintf('N_{surf}=%d / N_{tot}=%d', nR, Np)
    }, 'FontSize', 7);
end
sgtitle('Surface Diffuse Reflectance (log-scale)', 'FontSize', 14);

%% 2) Transmitted maps
figure('Color','w','Position',[100 100 1800 900]);
for k = 1:nFiles
    subplot(rows, cols, k);
    Tmap = Trans_all{k}.map;
    contourf(xEdges(1:end-1), yEdges(1:end-1), Tmap, 30, 'LineColor','none');
    axis equal tight; set(gca,'YDir','normal'); colorbar;
    clim([-1 1]);
    xlabel('x [cm]'); ylabel('y [cm]');

    params = parseFilename(files(k).name);
    nT = Trans_all{k}.nPhot; Np = Trans_all{k}.Nphot;

    title({
        sprintf('us=[%.3g %.3g %.3g]', params.us1, params.us2, params.us3)
        sprintf('ua5=%.3g us5=%.3g n5=%.2f', params.ua5, params.us5, params.n5)
        sprintf('N_{trans}=%d / N_{tot}=%d', nT, Np)
    }, 'FontSize', 7);
end
sgtitle('Transmitted Photons (bottom+sides, log-scale)', 'FontSize', 14);

%% 3) Trapped-in-ink maps
figure('Color','w','Position',[100 100 1800 900]);
for k = 1:nFiles
    subplot(rows, cols, k);
    Imap = Ink_all{k}.map;
    contourf(xEdges(1:end-1), yEdges(1:end-1), Imap, 30, 'LineColor','none');
    axis equal tight; set(gca,'YDir','normal'); colorbar;
    colormap turbo;
    clim([-20 1]);
    xlabel('x [cm]'); ylabel('y [cm]');
    xlim([-5 5]); ylim([-5 5]);

    params = parseFilename(files(k).name);
    nI = Ink_all{k}.nPhot; Np = Ink_all{k}.Nphot;

    title({
        sprintf('us=[%.3g %.3g %.3g]', params.us1, params.us2, params.us3)
        sprintf('ua5=%.3g us5=%.3g n5=%.2f', params.ua5, params.us5, params.n5)
        sprintf('N_{ink}=%d / N_{tot}=%d', nI, Np)
    }, 'FontSize', 7);
end
sgtitle('Photons Trapped in Ink (log-scale)', 'FontSize', 14);

%% gradient model diffuse reflectance at x axis, y=0

ks = [1 6 11];

customNames = {
    'No dye'
    '1 layer dye'
    '2 layer dye'
};

figure('Color','w'); hold on;

xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
[X,Y] = meshgrid(xCenters, yCenters);
Rrad  = sqrt(X.^2 + Y.^2);

dr = mean(diff(xCenters));
rBins = 0:dr:max(xCenters);
rCenters = (rBins(1:end-1) + rBins(2:end)) / 2;

for ii = 1:numel(ks)
    kk = ks(ii);

    RRlog = Diff_refl{kk}.map;   % log10(R)
    RRlin = 10.^RRlog;           % linear reflectance

    R_mean = nan(size(rCenters));
    for i = 1:numel(rCenters)
        mask = (Rrad >= rBins(i)) & (Rrad < rBins(i+1));
        vals = RRlin(mask);
        if ~isempty(vals)
            R_mean(i) = mean(vals,'omitnan');
        end
    end

    plot(rCenters, R_mean, ...
        'LineWidth', 1.5, ...
        'DisplayName', customNames{ii});
end

set(gca,'YScale','log');
xlabel('Distance r (cm)');
ylabel('Reflectance R (a.u.)');
title('MC radial reflectance (linear average, then log)');
legend('Location','best');
grid on;



%% No dye / dye 1 / dye 2 layers
figure; clf; set(gcf,'Color','w','Position',[100 100 1800 900]);

ks = [1 6 11];
customNames = {
    'No dye'
    '1 layer dye'
    '2 layer dye'
};

rows = 1; cols = numel(ks);

for i = 1:numel(ks)
    k = ks(i);
    subplot(rows, cols, i);

    RRlog = Diff_refl{k}.map;     % log10(R)

    % --- choose what you want to display ---
    RR = 10.^RRlog;               % linear R for plotting
    % RR = RRlog;                 % uncomment this instead if you want log10(R)

    contourf(xEdges(1:end-1), yEdges(1:end-1), RR', 30, 'LineColor','none');
    axis equal tight;
    set(gca,'YDir','normal');

    cb = colorbar;
    colormap hot;

    % colorbar label consistent with what you plot
    ylabel(cb, 'Diffuse reflectance R (a.u.)');      % for linear RR
    % ylabel(cb, 'log_{10}(R)');                     % if you plot RRlog instead

    % --- set color limits consistent with what you plot ---
    clim([0 300]);     % keep yours (only makes sense for linear RR)
    % clim([-5 0]);    % example if plotting RRlog

    xlabel('x (cm)');
    ylabel('y (cm)');
    xlim([-3.5 3.5]);
    ylim([-3.5 3.5]);

    nR = Diff_refl{k}.nPhot;
    Np = Diff_refl{k}.Nphot;

    title({
        customNames{i}
        sprintf('N_{surf} = %d / N_{tot} = %d', nR, Np)
    },'Fontname', 'Arial', 'FontSize', 10);
end

sgtitle('Surface Diffuse Reflectance', 'FontSize', 14);



%% plot mean pathlengths in each layers 
L0 = LayerPath{1}.Li;
L1 = LayerPath{6}.Li;
L2 = LayerPath{11}.Li;

figure;

B = [L0(:), L1(:), L2(:)];   % layers × conditions
bar(B);

legend({'No dye','1-layer dye','2-layer dye'}, 'Location','best');
xlabel('Layer');
ylabel('Mean photon pathlength [cm]');

set(gca,'XTickLabel',{'L1','L2','L3','Bulk'});
title('Photon pathlength redistribution due to dye');

grid on;


%% 3) Trapped-in-ink PHOTON COUNT maps (files 5, 10, 15 only)

filesToPlot = [5 10 15];
figure; clf; set(gcf,'Color','w','Position',[100 100 1800 900]);

customNames = {
    'No dye, \mu_a = 100 cm^{-1}'
    '1 layer dye, \mu_a = 100 cm^{-1}'
    '2 layer dye, \mu_a = 100 cm^{-1}'
};

for i = 1:numel(filesToPlot)

    k = filesToPlot(i);          % actual file index
    subplot(1,3,i);

    fileLabel = customNames{i};  % ✅ FIXED

    % --- Load file ---
    S = load(fullfile(files(k).folder, files(k).name));
    scatter_stat = S.scatter_stat;
    t_model      = S.t_model;

    % --- Photon exit positions ---
    x_ot = scatter_stat(:,4);
    y_ot = scatter_stat(:,5);
    z_ot = scatter_stat(:,6);

    % --- Ink geometry (same as MC) ---
    ink_mask = (x_ot >= -7.5) & (x_ot <= 7.5) & ...
               ((y_ot - 0).^2 + (z_ot - 1.0).^2 <= 0.3^2);

    % --- Exclude escaped photons ---
    idx_surface = (z_ot <= zMin);
    idx_trans   = (z_ot >= zMax) | ...
                  (x_ot <= xMin) | (x_ot >= xMax) | ...
                  (y_ot <= yMin) | (y_ot >= yMax);

    ink_trapped = ink_mask & ~idx_surface & ~idx_trans;

    % --- Coordinates of trapped photons ---
    x_ink = x_ot(ink_trapped);
    y_ink = y_ot(ink_trapped);

    % --- Photon COUNT map ---
    CountMap = zeros(Nx, Ny);
    [~,~,~,ix,iy] = histcounts2(x_ink, y_ink, xEdges, yEdges);

    for p = 1:numel(ix)
        if ix(p) > 0 && iy(p) > 0
            CountMap(ix(p), iy(p)) = CountMap(ix(p), iy(p)) + 1;
        end
    end

    % --- Plot ---
    contourf(xEdges(1:end-1), yEdges(1:end-1), CountMap', ...
             30, 'LineColor','none');
    axis equal tight;
    set(gca,'YDir','normal');
    colormap turbo;
    colorbar;
    caxis([0 2000]);

    xlabel('x [cm]');
    ylabel('y [cm]');
    xlim([-2 2]);
    ylim([-2 2]);

    % --- Title info ---
    nI = sum(ink_trapped);
    Np = size(scatter_stat,1);

    title({
        fileLabel
        sprintf('N_{ink}/N_{tot} = %.2d', nI/ Np)
    },'FontName','Arial','FontSize', 10);

end

%sgtitle('Photon count trapped inside ink region', 'FontSize', 14);


%% 4) Banana / voxelized fluence maps (3x5) — NEW (robust tiledlayout)
if plotBanana && any(~cellfun(@isempty, BananaSlices))

    % shared color limits for comparability
    allVals = [];
    for k = 1:nFiles
        if ~isempty(BananaSlices{k})
            allVals = [allVals; BananaSlices{k}(:)]; %#ok<AGROW>
        end
    end
    climB  = [prctile(allVals,5), prctile(allVals,99.5)];
    levels = linspace(climB(1), climB(2), bananaLevels);

    % ---- FIX: make sure layout has enough tiles ----
    rowsB = nRows;   % should be 3
    colsB = nCols;   % should be 5
    if rowsB*colsB < nFiles
        % fallback if nFiles > 15 or rows/cols got changed
        colsB = ceil(sqrt(nFiles));
        rowsB = ceil(nFiles/colsB);
    end

    figure('Color','w','Position',[100 100 1800 900]);
    t = tiledlayout(rowsB, colsB, 'Padding','compact', 'TileSpacing','compact');

    for k = 1:nFiles
        nexttile(t);   % ---- FIX: explicitly use this layout ----

        if isempty(BananaSlices{k})
            text(0.1,0.5,'No banana in file','Units','normalized');
            axis off; 
            continue;
        end

        contourf(BananaSlices{k}, levels, 'LineColor','none'); %
        axis image tight;
        xlim([50 250]);ylim([50 100]);
        caxis([-15 6]);
        colorbar off;

        params = parseFilename(files(k).name);
        title({
            sprintf('us=[%.3g %.3g %.3g]', params.us1, params.us2, params.us3)
            sprintf('ua5=%.3g us5=%.3g n5=%.2f', params.ua5, params.us5, params.n5)
        }, 'FontSize', 7);
    end

    cb = colorbar; 
    %colormap perula;
    cb.Layout.Tile = 'east';
    sgtitle(sprintf('Voxelized banana: log_{10}(%s) (XZ @ y≈0, rotated 180°)', bananaUse), 'FontSize', 14);
end



%% =========================================================
%  XZ banana profiles for cases 1, 6, 11 (center Y)
% =========================================================



ks = [1 6 11];   % no dye, 1-layer dye, 2-layer dye
figure('Color','w','Position',[100 100 1400 400]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

customNames = {
    'No dye'
    '1 layer dye'
    '2 layer dye'
};

for ii = 1:numel(ks)
    k = ks(ii);

    % ---- load banana data ----
    S = load(fullfile(files(k).folder, files(k).name), ...
             'B_wds','t_model');

    B = S.B_wds;                 % nx x ny x nz

    % ---- grid info ----
    nx = size(B,1);
    ny = size(B,2);
    nz = size(B,3);

    Lx = S.t_model.G.Lx;
    Lz = S.t_model.G.Lz;

    xCenters = ((1:nx)-0.5)*(Lx/nx) - Lx/2;
    zCenters = ((1:nz)-0.5)*(Lz/nz);

    % ---- center Y slice ----
    iy0 = round(ny/2);

    % ---- XZ slice using squeeze ----
    Bxz = squeeze(B(:,iy0,:))';   % nz x nx

    % ---- log + rotate 180° ----
    %Bxz = Bxz / max(Bxz(:));
    Bxz = flipud(fliplr(log10(Bxz + eps)));

    % ---- plot ----
   nexttile;
contourf(xCenters, zCenters, Bxz, 30, 'LineColor','none');
axis image;
set(gca,'YDir','normal');
colormap turbo;

title(customNames{ii}, 'Interpreter','tex');

xlim([-3 3]);
ylim([2.5 5]);
xlabel('x (cm)');
ylabel('z (cm)');

end

sgtitle('Banana profiles', 'FontSize', 14);
set(gca,'FontName','Arial','FontSize', 10);


%% Differential XZ banana: ink absorption effect (only files 5,10,15)
dy = 2;
filesToPlot = [5 10 15];

customNames = {
    'No dye, \mu_a = 100 cm^{-1}'
    '1 layer dye, \mu_a = 100 cm^{-1}'
    '2 layer dye, \mu_a = 100 cm^{-1}'
};

% --- physical dimensions (cm) ---
Lx = 15;   % cm  (so x goes -7.5 to +7.5)
Lz = 5;    % cm  (depth 0 to 5)

figure('Color','w','Position',[100 100 1800 700]);
t = tiledlayout(1, numel(filesToPlot), 'Padding','compact','TileSpacing','compact');

for i = 1:numel(filesToPlot)
    k = filesToPlot(i);

    % --- row-based reference ---
    row   = ceil(k / nCols);
    k_ref = (row-1)*nCols + 1;

    Sref = load(fullfile(files(k_ref).folder, files(k_ref).name), 'B_wds');
    Bref = Sref.B_wds;

    S = load(fullfile(files(k).folder, files(k).name), 'B_wds');
    B = S.B_wds;

    % center-Y slab
    ny  = size(B,2);
    iy0 = round(ny/2);
    ys  = max(1,iy0-dy) : min(ny,iy0+dy);

    % slab-averaged XZ
    Bref_xz = squeeze(mean(Bref(:,ys,:),2));
    B_xz    = squeeze(mean(B(:,ys,:),2));

    % log-difference
    Delta = log10(B_xz + eps) - log10(Bref_xz + eps);

    % orientation (keep your preferred orientation)
    Delta = rot90(Delta, -1);
    Delta = flipud(Delta);

    % --- centered physical axes ---
    Nz = size(Delta,1);   % depth pixels
    Nx = size(Delta,2);   % lateral pixels

    x = linspace(-Lx/2, Lx/2, Nx);   % cm, source at x=0
    z = linspace(0,    Lz,    Nz);   % cm, surface at z=0

    nexttile(t);
    contourf(x, z, Delta, 20, 'LineColor','none');
    axis image;
    set(gca,'YDir','normal');     % depth increases downward
    caxis([-1 0]);

    % make axes uniform across tiles
    xlim([-Lx/2, Lx/2]);
    ylim([0, Lz]);

    xlabel('Lateral distance x (cm)');
    if i == 1
        ylabel('Depth z (cm)');
    end

    title(customNames{i}, 'Interpreter','tex');
end

cb = colorbar;
cb.Label.String = '\Delta log_{10}(Intensity)';
cb.Layout.Tile = 'east';
set(gca,'FontName','Arial','FontSize', 10);
sgtitle('Photon intensity for ink–dye effect');

%% =========================================================
%  Depth-wise fluence F(z) for cases 1, 6, 11
% =========================================================

ks = [1 6 11];   % no dye, 1-layer dye, 2-layer dye
figure; hold on; set(gcf,'Color','w');

for ii = 1:numel(ks)
    k = ks(ii);

    % ---- load banana data ----
    S = load(fullfile(files(k).folder, files(k).name), 'B_wds', 't_model');

    B = S.B_wds;                     % nx x ny x nz
    nz = size(B,3);

    % ---- depth axis (cm) ----
    Lz = S.t_model.G.Lz;
    zCenters = ((1:nz) - 0.5) * (Lz / nz);

    % ---- depth-integrated fluence ----
    Fz = squeeze(sum(sum(B,1),2));    % sum over x,y
    %Fz = Fz / max(Fz);                % normalize

    plot(zCenters, Fz, 'LineWidth', 2, ...
        'DisplayName', sprintf('Case %d', k));
end

xlabel('Depth z (cm)');
ylabel('Normalized fluence F(z)');
title('Depth-wise fluence profiles (banana integral)');
legend('Location','best');
grid on;


%% contrast vs dye depth for an ink strength (shared axes)

figure('Color','w','Position',[100 100 1600 900]);

t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

inkCases = {
    [2 7 12],  [1 6 11],  'Ink \mu_a = 1 cm^{-1}'
    [3 8 13],  [1 6 11],  'Ink \mu_a = 7 cm^{-1}'
    [4 9 14],  [1 6 11],  'Ink \mu_a = 30 cm^{-1}'
    [5 10 15], [1 6 11],  'Ink \mu_a = 100 cm^{-1}'
};

colors = lines(3);
labels = {'No dye','1 layer dye','2 layer dye'};

ax = gobjects(4,1);   % store axes handles

for s = 1:4
    ax(s) = nexttile; hold on;

    sigIdx  = inkCases{s,1};
    baseIdx = inkCases{s,2};

    for i = 1:3
        [r, Rs] = radialAverage(Diff_refl{sigIdx(i)}.map,  xEdges, yEdges);
        [~, Rb] = radialAverage(Diff_refl{baseIdx(i)}.map, xEdges, yEdges);

        Contrast = Rs - Rb;

        plot(r, Contrast, 'LineWidth', 1, ...
            'Color', colors(i,:), ...
            'DisplayName', labels{i});
    end

    title(inkCases{s,3});
    yline(0,'k--');
    grid off;

    xlim([0 4]);
    ylim([-0.6 0.2]);
end

% ---- Link all axes ----
linkaxes(ax,'xy');

% ---- Shared labels ----
xlabel(t,'Distance (cm)','FontSize',12);
ylabel(t,'Signal contrast (a.u.)','FontSize',12);

legend(ax(1),'Location','best');
set(gca,'FontName','Arial','FontSize', 10);
%sgtitle('Effect of dye depth on ink contrast (radial average)', ...
%        'FontSize',14);


%%
function [rCenters, Ravg] = radialAverage(RR, xEdges, yEdges)
% RR : Nx × Ny map (log10 reflectance)
% Returns azimuthally averaged R(r)

xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;

[X,Y] = meshgrid(xCenters, yCenters);
Rmap  = RR';                       % match X,Y orientation
Rrad  = sqrt(X.^2 + Y.^2);

dr = mean(diff(xCenters));
rBins = 0:dr:max(xCenters);
rCenters = (rBins(1:end-1) + rBins(2:end)) / 2;

Ravg = nan(size(rCenters));

for i = 1:numel(rCenters)
    mask = (Rrad >= rBins(i)) & (Rrad < rBins(i+1));
    vals = Rmap(mask);
    if ~isempty(vals)
        Ravg(i) = mean(vals,'omitnan');
    end
end
Ravg = Ravg;
%Ravg = smoothdata(Ravg,'sgolay',1);
end

%% ===================== END MAIN ==========================

%% =========================================================
%  Filename parser (UPDATED for "p" decimals)
% =========================================================
function params = parseFilename(fname)
tokens = regexp(fname, ...
    'us1_(.*?)_us2_(.*?)_us3_(.*?)_ua5_(.*?)_us5_(.*?)_n5_(.*?).mat', ...
    'tokens');
tokens = tokens{1};

toNum = @(s) str2double(strrep(s,'p','.'));  % 0p8 -> 0.8

params.us1 = toNum(tokens{1});
params.us2 = toNum(tokens{2});
params.us3 = toNum(tokens{3});
params.ua5 = toNum(tokens{4});
params.us5 = toNum(tokens{5});
params.n5  = toNum(tokens{6});
end
