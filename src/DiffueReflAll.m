%% =========================================================
%  Diffuse Reflectance Surface Maps for ALL Conditions
%  Uses scatter_stat (final photon stats)
%  Author: Arnab Paul
% ==========================================================
clear; clc;

% -------------------- USER SETTINGS ----------------------
dataFolder = pwd;      % folder containing all .mat files
tol = 0.01;            % surface tolerance (cm) ~0.1 mm
zSurface = 0;          % detector surface at z = 0
useLog = true;         % log-scale plotting
% ----------------------------------------------------------

%% --------- BUILD ORDERED FILE LIST (3x5 grid) ----------
allFiles   = dir(fullfile(dataFolder, '*.mat'));
nAll       = numel(allFiles);

us1 = zeros(nAll,1);
ua3 = zeros(nAll,1);

for k = 1:nAll
    % sim_info must contain us1 and ua3 for this to work
    S = load(fullfile(allFiles(k).folder, allFiles(k).name), 'sim_info');
    us1(k) = S.sim_info.us1;   % dye condition: 30=No, 10=1L, 2=2L
    ua3(k) = S.sim_info.ua3;   % ink condition: 0.05,1,7,30,100
end

% Desired row/column values (fix order here)
rowUs1 = [30, 10, 2];              % rows: No dye, 1L dye, 2L dye
colUa3 = [0.05, 1, 7, 30, 100];    % cols: No ink, 1, 7, 30, 100

nRows = numel(rowUs1);
nCols = numel(colUa3);

% Build ordered files struct: 3×5 in the desired serial order
ordered = repmat(allFiles(1), nRows*nCols, 1);
idx = 0;
for i = 1:nRows
    for j = 1:nCols
        idx = idx + 1;
        mask = (us1 == rowUs1(i)) & (ua3 == colUa3(j));
        fidx = find(mask, 1);
        if isempty(fidx)
            error('No file for us1=%g, ua3=%g', rowUs1(i), colUa3(j));
        end
        ordered(idx) = allFiles(fidx);
    end
end

files  = ordered;
nFiles = numel(files);
rows   = nRows;
cols   = nCols;

figure('Color','w','Position',[100 100 1800 900]);
R_all = cell(nFiles,1);      % or zeros(Ny,Nx,nFiles) if Ny,Nx fixed
%% -------------------- MAIN LOOP --------------------------
for k = 1:nFiles
    % Load file
    S = load(fullfile(files(k).folder, files(k).name));
    scatter_stat = S.scatter_stat;
    t_model = S.t_model;

    % ----- Extract final photon data -----
    x_ot = scatter_stat(:,4);
    y_ot = scatter_stat(:,5);
    z_ot = scatter_stat(:,6);
    w_ot = scatter_stat(:,8);

    % ----- Select surface photons -----
    idx_surface = abs(z_ot - zSurface) <= tol;
    x_surf = x_ot(idx_surface);
    y_surf = y_ot(idx_surface);
    w_surf = w_ot(idx_surface);

    % ----- Grid info -----
    Lx = t_model.G.Lx;
    Ly = t_model.G.Ly;
    Nx = t_model.G.nx;
    Ny = t_model.G.ny;
    xEdges = linspace(-Lx/2, Lx/2, Nx+1);
    yEdges = linspace(-Ly/2, Ly/2, Ny+1);

    % ----- Bin photons -----
    xi = discretize(x_surf, xEdges);
    yi = discretize(y_surf, yEdges);
    valid = ~isnan(xi) & ~isnan(yi);
    xi = xi(valid);
    yi = yi(valid);
    w  = w_surf(valid);

    % ----- Reflectance map -----
    R = accumarray([yi, xi], w, [Ny, Nx], @sum, 0);
    if useLog
        R(R <= 0) = NaN;
        R = log10(R);
    end
    R_all{k} = R;            % <--- store reflectance map
    %k=(row−1)×5+col; so, for example, row 3 col 1 (Dye2, no signal) is k = (3-1)*5 + 1 = 11, so R_all{11}
    %% ----- Plot -----
    subplot(rows, cols, k)
    contourf(xEdges(1:end-1), yEdges(1:end-1), R, 30,'LineColor','none'); % add ,'LineColor','none' if desired
    axis equal tight
    set(gca,'YDir','normal')
    %colormap hot
    colorbar
    clim([-5 5]);
    xlabel('x [cm]')
    ylabel('y [cm]')
   %% plot Diffuse reflectance improvements after dye
    R_nd   = R_all{1};
    R_d1   = R_all{6};
    R_d2   = R_all{11};
    R_no_lin  = 10.^R_nd;     % back to linear
    R_d1_lin  = 10.^R_d1;
    R_d2_lin  = 10.^R_d2;
    

    diff_d1_nd = R_d1_lin - R_no_lin;
    diff_d2_nd = R_d2_lin - R_no_lin;
    figure;
    subplot(1,2,1); surf(diff_d1_nd); axis image; colorbar;clim([0 5]);
    title('Dye1 - No dye (Nosignal)');
    
    subplot(1,2,2); surf(diff_d2_nd); axis image; colorbar;clim([0 5]);
    title('Dye2 - No dye (Nosignal)');
    
    colormap hot;


    %% ----- Title from filename -----
    params = parseFilename(files(k).name);
    title({
        sprintf('L1: ua=%.3g  us=%.3g  n=%.2f', params.ua1, params.us1, params.n1)
        sprintf('L2: ua=%.3g  us=%.3g  n=%.2f', params.ua2, params.us2, params.n2)
        sprintf('L3: ua=%.3g  us=%.3g  n=%.2f', params.ua3, params.us3, params.n3)
        }, 'FontSize', 7);
end

sgtitle('Surface Diffuse Reflectance (log-scale)', 'FontSize', 14);

%% ===================== END MAIN ==========================

%% =========================================================
%  Filename parser (DO NOT MODIFY)
% =========================================================
function params = parseFilename(fname)
tokens = regexp(fname, ...
    'ua1_(.*?)_us1_(.*?)_n1_(.*?)_ua2_(.*?)_us2_(.*?)_n2_(.*?)_ua3_(.*?)_us3_(.*?)_n3_(.*?).mat', ...
    'tokens');
tokens = tokens{1};
params.ua1 = str2double(tokens{1});
params.us1 = str2double(tokens{2});
params.n1  = str2double(tokens{3});
params.ua2 = str2double(tokens{4});
params.us2 = str2double(tokens{5});
params.n2  = str2double(tokens{6});
params.ua3 = str2double(tokens{7});
params.us3 = str2double(tokens{8});
params.n3  = str2double(tokens{9});
end
