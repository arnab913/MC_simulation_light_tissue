%% =========================================================
%  Diffuse Reflectance / Transmission / Ink Maps for ALL Conditions
%  Uses scatter_stat (final photon stats)
% ==========================================================
clear; clc;

% -------------------- USER SETTINGS ----------------------
%dataFolder = pwd;      % folder containing all .mat files
dataFolder = ("E:\MC SImulation data\Data_Gradient_tissue_123cmlayers");
xMin = -7.5; xMax = 7.5;
yMin = -7.5; yMax = 7.5;
zMin = 0;    zMax = 5;      % top z = 0, bottom z = 5 (adjust if needed)
useLog = true;              % log-scale maps
% ----------------------------------------------------------

%% --------- BUILD ORDERED FILE LIST (3x5 grid) ----------
allFiles   = dir(fullfile(dataFolder, '*.mat'));
nAll       = numel(allFiles);

us1 = zeros(nAll,1);
ua5 = zeros(nAll,1);
for k = 1:nAll
    S = load(fullfile(allFiles(k).folder, allFiles(k).name), 'sim_info');
    us1(k) = S.sim_info.us1;
    ua5(k) = S.sim_info.ua5;
end

rowUs1 = [30, 5, 2];              % rows: No dye, 1L dye, 2L dye
colUa5 = [0.05, 1, 7, 30, 150];    % cols: No ink, 1, 7, 30, 100
nRows = numel(rowUs1);
nCols = numel(colUa5);

ordered = repmat(allFiles(1), nRows*nCols, 1);
idx = 0;
for i = 1:nRows
    for j = 1:nCols
        idx = idx + 1;
        mask = (us1 == rowUs1(i)) & (ua5 == colUa5(j));
        fidx = find(mask, 1);
        if isempty(fidx)
            error('No file for us1=%g, ua3=%g', rowUs1(i), colUa5(j));
        end
        ordered(idx) = allFiles(fidx);
    end
end

files  = ordered;
nFiles = numel(files);

Diff_refl = cell(nFiles,1);   % each cell: struct with map + counts
Trans_all = cell(nFiles,1);
Ink_all   = cell(nFiles,1);

%% -------------------- MAIN LOOP --------------------------
for k = 1:nFiles

    % ----- Load file -----
    S = load(fullfile(files(k).folder, files(k).name));
    scatter_stat = S.scatter_stat;
    t_model      = S.t_model;

    % ----- Extract final photon data -----
    x_ot = scatter_stat(:,4);
    y_ot = scatter_stat(:,5);
    z_ot = scatter_stat(:,6);
    w_ot = scatter_stat(:,8);

    % total launched photons (assume one row per photon)
    Nphot = size(scatter_stat,1);

    % ----- Grid info from model -----
    Lx = t_model.G.Lx;
    Ly = t_model.G.Ly;
    Lz = t_model.G.Lz;
    Nx = t_model.G.nx;
    Ny = t_model.G.ny;
    %Nz = t_model.G.nz;   % not used here

    xEdges = linspace(-Lx/2, Lx/2, Nx+1);
    yEdges = linspace(-Ly/2, Ly/2, Ny+1);
    %zEdges = linspace( 0, Lz, Nz+1);  % not used for 2D maps

    % --------------------------------------------------
    % 1) Surface photons (diffuse reflectance, top at zMin)
    % --------------------------------------------------
    idx_surface = (z_ot <= zMin);

    x_surf = x_ot(idx_surface);
    y_surf = y_ot(idx_surface);
    w_surf = w_ot(idx_surface);

    n_reflect = sum(idx_surface);        % count

    % Bin reflectance on x-y plane
    [~,~,~,idxR,idyR] = histcounts2(x_surf, y_surf, xEdges, yEdges);
    RR = zeros(Nx,Ny);
    for i = 1:numel(idxR)
        if idxR(i) > 0 && idyR(i) > 0
            RR(idxR(i), idyR(i)) = RR(idxR(i), idyR(i)) + w_surf(i);
        end
    end
    if useLog
        RR(RR <= 0) = NaN;
        RR = log10(RR);
    end

    % store as struct in cell
    Diff_refl{k} = struct( ...
        'map',   RR, ...
        'nPhot', n_reflect, ...
        'Nphot', Nphot);

    % --------------------------------------------------
    % 2) Transmitted photons (bottom + sides, no top)
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

    n_trans = sum(idx_transmitted);      % count

    [~,~,~,idxT,idyT] = histcounts2(x_trans, y_trans, xEdges, yEdges);
    Trans_map = zeros(Nx,Ny);
    for i = 1:numel(idxT)
        if idxT(i) > 0 && idyT(i) > 0
            Trans_map(idxT(i), idyT(i)) = Trans_map(idxT(i), idyT(i)) + w_trans(i);
        end
    end
    if useLog
        Trans_map(Trans_map <= 0) = NaN;
        Trans_map = log10(Trans_map);
    end

    Trans_all{k} = struct( ...
        'map',   Trans_map, ...
        'nPhot', n_trans, ...
        'Nphot', Nphot);

    % --------------------------------------------------
    % 3) Photons trapped in the ink layer (not escaped)
    % --------------------------------------------------
    ink_mask = (x_ot >= -7.5) & (x_ot <= 7.5) & ...
               ((y_ot - 0).^2 + (z_ot - 1.0).^2 <= 0.3^2);
    ink_trapped_mask = ink_mask & ~idx_surface & ~idx_transmitted;

    x_ink = x_ot(ink_trapped_mask);
    y_ink = y_ot(ink_trapped_mask);
    w_ink = w_ot(ink_trapped_mask);

    n_ink = sum(ink_trapped_mask);       % count

    [~,~,~,idxI,idyI] = histcounts2(x_ink, y_ink, xEdges, yEdges);
    Ink_map = zeros(Nx,Ny);
    for i = 1:numel(idxI)
        if idxI(i) > 0 && idyI(i) > 0
            Ink_map(idxI(i), idyI(i)) = Ink_map(idxI(i), idyI(i)) + w_ink(i);
        end
    end
    if useLog
        Ink_map(Ink_map <= 0) = NaN;
        Ink_map = log10(Ink_map);
    end

    Ink_all{k} = struct( ...
        'map',   Ink_map, ...
        'nPhot', n_ink, ...
        'Nphot', Nphot);

    % --------------------------------------------------
    % 4) Remaining photons (inside tissue, not in any category)
    %     (optional: counts only)
    % --------------------------------------------------
    idx_classified = idx_surface | idx_transmitted | ink_trapped_mask;
    idx_tissue     = ~idx_classified;
    n_tissue       = sum(idx_tissue);    %#ok<NASGU>  % if you want to inspect later




%% Mean pathlength per layer
% --------------------------------------------------
% 5) Mean photon pathlength per layer (Li estimation)
% --------------------------------------------------

% Final depth of photons
z_out = scatter_stat(:,6);
s_tot = scatter_stat(:,7);   % total pathlength per photon

% Define layer boundaries (cm)
zL = [0.0, 0.3, 0.4, 0.5, inf];   % 3 dye layers + bulk

nLayers = numel(zL)-1;
Li = zeros(nLayers,1);          % accumulated pathlength
Ni = zeros(nLayers,1);          % photon counts

for p = 1:Nphot

    zf = z_out(p);
    s  = s_tot(p);

    % Skip photons that never entered tissue (rare)
    if s <= 0
        continue
    end

    % Vertical distance traversed in each layer
    dz = zeros(nLayers,1);

    for j = 1:nLayers
        z1 = zL(j);
        z2 = zL(j+1);

        % overlap of photon path with layer
        dz(j) = max(0, min(zf, z2) - z1);
    end

    dz_tot = sum(dz);

    if dz_tot > 0
        % distribute total pathlength proportionally
        Li = Li + s * (dz / dz_tot);
        Ni = Ni + (dz > 0);
    end
end

% Mean pathlength per layer
Li_mean = Li ./ max(Ni,1);

% Store
LayerPath{k} = struct( ...
    'Li', Li_mean, ...
    'zBounds', zL );

end




%% ================== PLOTTING ============================
% grid size for subplots
rows = nRows;
cols = nCols;

% Rebuild grid (same as in loop) from last t_model
Lx = t_model.G.Lx;
Ly = t_model.G.Ly;
Nx = t_model.G.nx;
Ny = t_model.G.ny;
xEdges = linspace(-Lx/2, Lx/2, Nx+1);
yEdges = linspace(-Ly/2, Ly/2, Ny+1);

%% 1) Diffuse Reflectance maps
figure; clf; set(gcf,'Color','w','Position',[100 100 1800 900]);
for k = 1:nFiles
    subplot(rows, cols, k);
    RR = Diff_refl{k}.map;
    contourf(xEdges(1:end-1), yEdges(1:end-1), RR, 30, 'LineColor','none');
    axis equal tight;
    set(gca,'YDir','normal');
    colorbar;
    clim([-5 5]);
    xlabel('x [cm]');
    ylabel('y [cm]');

       params = parseFilename(files(k).name);
    nR = Diff_refl{k}.nPhot;
    Np = Diff_refl{k}.Nphot;

    title({
    sprintf('us1=%.3g  us2=%.3g  us3=%.3g', params.us1, params.us2, params.us3)
    sprintf('ua5=%.3g  us5=%.3g  n5=%.2f', params.ua5, params.us5, params.n5)
    sprintf('N_{surf} = %d / N_{tot} = %d', nR, Np)
    }, 'FontSize', 7);

end
sgtitle('Surface Diffuse Reflectance (log-scale)', 'FontSize', 14);

%% 2) Transmitted maps (bottom + sides)
figure; clf; set(gcf,'Color','w','Position',[100 100 1800 900]);
for k = 1:nFiles
    subplot(rows, cols, k);
    Tmap = Trans_all{k}.map;
    contourf(xEdges(1:end-1), yEdges(1:end-1), Tmap, 30, 'LineColor','none');
    axis equal tight;
    set(gca,'YDir','normal');
    colorbar;
    clim([-1 1]);
    xlabel('x [cm]');
    ylabel('y [cm]');

       params = parseFilename(files(k).name);
    nT = Trans_all{k}.nPhot;
    Np = Trans_all{k}.Nphot;

    title({
    sprintf('us1=%.3g  us2=%.3g  us3=%.3g', params.us1, params.us2, params.us3)
    sprintf('ua5=%.3g  us5=%.3g  n5=%.2f', params.ua5, params.us5, params.n5)
    sprintf('N_{surf} = %d / N_{tot} = %d', nT, Np)
    }, 'FontSize', 7);

end
sgtitle('Transmitted Photons (bottom + sides, log-scale)', 'FontSize', 14);

%% 3) Trapped-in-ink maps
figure; clf; set(gcf,'Color','w','Position',[100 100 1800 900]);
for k = 1:nFiles
    subplot(rows, cols, k);
    Imap = Ink_all{k}.map;
    contourf(xEdges(1:end-1), yEdges(1:end-1), Imap, 30, 'LineColor','none');
    axis equal tight;
    set(gca,'YDir','normal');
    colorbar;
    colormap turbo;
    clim([-20 1]);
    xlabel('x [cm]'); xlim([-5 5]);
    ylabel('y [cm]'); ylim([-5 5]);

        params = parseFilename(files(k).name);
    nI = Ink_all{k}.nPhot;
    Np = Ink_all{k}.Nphot;

    title({
    sprintf('us1=%.3g  us2=%.3g  us3=%.3g', params.us1, params.us2, params.us3)
    sprintf('ua5=%.3g  us5=%.3g  n5=%.2f', params.ua5, params.us5, params.n5)
    sprintf('N_{surf} = %d / N_{tot} = %d', nI, Np)
    }, 'FontSize', 7);

end
sgtitle('Photons Trapped in Ink (log-scale)', 'FontSize', 14);

%% gradient model diffuse reflectance at x axis, y=0
ks = [1 6 11];
figure; hold on;

xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
[X,Y] = meshgrid(xCenters, yCenters);
Rrad  = sqrt(X.^2 + Y.^2);

dr = mean(diff(xCenters));
rBins = 0:dr:max(xCenters);
rCenters = (rBins(1:end-1) + rBins(2:end)) / 2;

for kk = ks
    RRlog = Diff_refl{kk}.map;     % log10(R)
    RRlin = 10.^RRlog;             % linear reflectance per pixel-bin

    R_mean = nan(size(rCenters));
    for i = 1:numel(rCenters)
        mask = (Rrad >= rBins(i)) & (Rrad < rBins(i+1));
        vals = RRlin(mask);
        if ~isempty(vals)
            R_mean(i) = mean(vals,'omitnan');    % mean linear reflectance
        end
    end

    plot(rCenters,R_mean, 'LineWidth', 1, 'DisplayName', sprintf('k=%d',kk));
   set(gca,'YScale','log')

end

xlabel('Distance r (cm)');
ylabel('log_{10} R (a.u.)');
title('MC radial reflectance (linear average, then log)');
legend show;
grid on;
%xlim([0 3.5]);




%% No dye dye 1 dye 2 layers
figure; clf; set(gcf,'Color','w','Position',[100 100 1800 900]);
ks   = [1 6 11];
rows = 1; cols = numel(ks);

for i = 1:numel(ks)
    k = ks(i);                   % actual case index
    subplot(rows, cols, i);

    RR = Diff_refl{k}.map;       % already log10(R)
    RR = 10.^RR;% for lineaar scaling, for log just comment this line

    contourf(xEdges(1:end-1), yEdges(1:end-1), RR', 30, 'LineColor','none');
    axis equal tight;
    set(gca,'YDir','normal');%,,'YScale','log'

    cb = colorbar;
    colormap hot;
    ylabel(cb, 'log_{10} diffuse reflectance');

    % Choose appropriate log range; adjust to your data
    % Example: reflectance between 10^-5 and 10^0 -> log10 in [-5, 0]
    clim([0 300]);    % or e.g. [-4 -1], depending on RR range

    xlabel('x [cm]');
    ylabel('y [cm]');
    xlim([-3.5 3.5]);
    ylim([-3.5 3.5]);

    nR = Diff_refl{k}.nPhot;
    Np = Diff_refl{k}.Nphot;
    title({
        sprintf('Case k = %d', k)
        sprintf('N_{surf} = %d / N_{tot} = %d', nR, Np)
    }, 'FontSize', 7);
end

sgtitle('Surface Diffuse Reflectance (log_{10} scale)', 'FontSize', 14);


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



%% ===================== END MAIN ==========================
% =========================================================
%  Filename parser (DO NOT MODIFY)
% =========================================================
function params = parseFilename(fname)
tokens = regexp(fname, ...
    'us1_(.*?)_us2_(.*?)_us3_(.*?)_ua5_(.*?)_us5_(.*?)_n5_(.*?).mat', ...
    'tokens');
tokens = tokens{1};
params.us1 = str2double(tokens{1});   % top tissue scattering
params.us2 = str2double(tokens{2});   % middle tissue scattering
params.us3 = str2double(tokens{3});   % bottom tissue scattering
params.ua5 = str2double(tokens{4});   % ink absorption in layer 5
params.us5 = str2double(tokens{5});   % ink scattering in layer 5
params.n5  = str2double(tokens{6});   % ink refractive index (or global n)
end

