%% =========================================================
%  Diffuse Reflectance Surface Maps for ALL Conditions
%  Uses scatter_stat (final photon stats)
%  Author: (you)
% ==========================================================

clear; clc;

% -------------------- USER SETTINGS ----------------------
dataFolder = pwd;      % folder containing all .mat files
tol = 0.01;            % surface tolerance (cm) ~0.1 mm
zSurface = 0;          % detector surface at z = 0
useLog = true;         % log-scale plotting
cols = 5;              % subplot columns (15 files → 3x5)
% ----------------------------------------------------------

files = dir(fullfile(dataFolder, '*.mat'));
nFiles = numel(files);

rows = ceil(nFiles / cols);

figure('Color','w','Position',[100 100 1800 900]);

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

    %% ----- Plot -----
    subplot(rows, cols, k)

    contourf(xEdges(1:end-1), yEdges(1:end-1), R, 30, 'LineColor','none');
    axis equal tight
    set(gca,'YDir','normal')
    colormap hot
    colorbar

    xlabel('x [cm]')
    ylabel('y [cm]')

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
