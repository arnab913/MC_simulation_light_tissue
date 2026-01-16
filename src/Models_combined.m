%% =========================================================
%  ALL MODELS IN ONE SCRIPT (MODEL-ONLY MC REFLECTANCE):
%   (1) MC radial reflectance recomputed for modeling:
%       R(x,y) = sum(w_surf)/(Nlaunch*dA)  (reflectance density)
%   (2) Farrell diffusion dipole (uses mu_s' = mu_s*(1-g))
%   (3) Schmitt layered model (REVERTED: mu_eff = sqrt(3*mua*mus) using saved us)
%  PLUS: overlay + subplots
%
%  REQUIRED IN WORKSPACE:
%    - files      : struct array from dir()/your ordered list; files(k).folder, files(k).name
%    - LayerPath  : cell array; LayerPath{k}.Li exists for k in ks
%
%  REQUIRED INSIDE EACH .mat FILE (files(k)):
%    - scatter_stat
%    - t_model
%    - sim_info
%
%  UNITS:
%    x,y,z in cm ; mua, mus in cm^-1
% =========================================================

% clearvars -except files LayerPath
% close all; clc;

%% -----------------------------
% USER SETTINGS
% -----------------------------
ks = [1 6 11];       % no dye, dye1, dye2 (file indices in your ordered list)

g = 0.90;            % anisotropy used in MC (if isotropic, set 0)

% MC surface criterion
zSurf = 0;           % cm (top surface)
zTol  = 1e-9;        % tolerance

% Farrell dipole settings
mua_farrell = 0.05;                   % [cm^-1]
n_farrell   = 1.5;                    % refractive index
r_farrell   = linspace(0.01, 6, 400); % [cm]

% You said simulation uses mu_s (NOT mu_s')
% So these are mu_s per layer [cm^-1], then convert to mu_s'
mus_layers = struct();
mus_layers.k1  = [30 30 30 30];
mus_layers.k6  = [5 10 20 30];
mus_layers.k11 = [2 5 10 30];

% Schmitt settings
labels_sch  = {'No dye','1 layer dye','2 layer dye'};
valid_r_min = 0.3;
valid_r_max = 4.0;

% Plot limits
xlim_plot = [0 6];
ylim_plot = [1e-5 10];

%% =========================================================
% (0) BASIC CHECKS
% =========================================================
if ~exist('files','var') || isempty(files)
    error('`files` not found. Provide your ordered file list (struct from dir).');
end
if ~exist('LayerPath','var') || isempty(LayerPath)
    error('`LayerPath` not found. Run your LayerPath-building code first.');
end
if max(ks) > numel(files)
    error('ks contains %d but only %d files available.', max(ks), numel(files));
end

%% =========================================================
% (1) MC radial reflectance FOR MODELS ONLY (recomputed)
%     R(x,y) = sum(w_surf) / (Nlaunch * dA)
% =========================================================
R_MC = [];        % will allocate after first case
rCenters = [];    % radial bins from model grid
xEdges = []; yEdges = [];  %#ok<NASGU> % kept for potential debugging

for a = 1:numel(ks)
    kk = ks(a);

    S = load(fullfile(files(kk).folder, files(kk).name), ...
        'scatter_stat','t_model');

    if ~isfield(S,'scatter_stat') || ~isfield(S,'t_model')
        error('scatter_stat or t_model missing inside %s', fullfile(files(kk).folder, files(kk).name));
    end

    scatter_stat = S.scatter_stat;
    t_model      = S.t_model;

    % columns (your convention)
    x = scatter_stat(:,4);
    y = scatter_stat(:,5);
    z = scatter_stat(:,6);
    w = scatter_stat(:,8);

    Nlaunch = size(scatter_stat,1);  % keep consistent with your previous pipeline

    % --- surface hits ---
    idxSurf = (z <= (zSurf + zTol));
    x_s = x(idxSurf);
    y_s = y(idxSurf);
    w_s = w(idxSurf);

    % --- grid from model ---
    Lx = t_model.G.Lx;
    Ly = t_model.G.Ly;
    Nx = t_model.G.nx;
    Ny = t_model.G.ny;

    xEdges_k = linspace(-Lx/2, Lx/2, Nx+1);
    yEdges_k = linspace(-Ly/2, Ly/2, Ny+1);

    dx = mean(diff(xEdges_k));
    dy = mean(diff(yEdges_k));
    dA = dx * dy;  % cm^2

    % --- accumulate weights to pixels ---
    [~,~,~,ix,iy] = histcounts2(x_s, y_s, xEdges_k, yEdges_k);

    RRw = zeros(Nx, Ny); % weight-sum per pixel
    for p = 1:numel(ix)
        if ix(p) > 0 && iy(p) > 0
            RRw(ix(p), iy(p)) = RRw(ix(p), iy(p)) + w_s(p);
        end
    end

    % --- diffusion-consistent reflectance density ---
    Rmap_lin = RRw ./ (Nlaunch * dA);  % [1/cm^2] (relative units)

    % --- radial average (linear) ---
    [rC, R_mean] = radialAverage_linMap(Rmap_lin, xEdges_k, yEdges_k);

    % initialize storage using first case bins
    if a == 1
        rCenters = rC(:).';                     % row
        R_MC = nan(numel(ks), numel(rCenters)); % rows for cases

        % keep edges for consistency check
        xEdges = xEdges_k; %#ok<NASGU>
        yEdges = yEdges_k; %#ok<NASGU>
    else
        % safety: ensure same binning across cases
        if numel(rC) ~= numel(rCenters) || any(abs(rC(:) - rCenters(:)) > 1e-12)
            error('Radial bins differ between cases. Ensure all cases share same grid Nx,Ny,Lx,Ly.');
        end
    end

    R_MC(a,:) = R_mean;

    fprintf('Model-MC k=%d: Nlaunch=%d, Nsurf=%d, dx=%.4g cm, dA=%.4g cm^2\n', ...
        kk, Nlaunch, sum(idxSurf), dx, dA);
end

%% =========================================================
% (2) Farrell diffusion dipole (compute & store)
% =========================================================
mu_sp_eff_all = zeros(size(ks));

for idx = 1:numel(ks)
    k = ks(idx);

    if numel(LayerPath) < k || isempty(LayerPath{k}) || ~isfield(LayerPath{k},'Li')
        error('LayerPath{%d}.Li not found.', k);
    end

    Li = LayerPath{k}.Li;
    if numel(Li) < 4
        error('LayerPath{%d}.Li must have at least 4 elements.', k);
    end
    Li = Li(1:4);
    Li = Li(:); % column

    switch k
        case 1,  mus_i = mus_layers.k1;
        case 6,  mus_i = mus_layers.k6;
        case 11, mus_i = mus_layers.k11;
        otherwise, error('Unexpected k=%d', k);
    end
    mus_i = mus_i(:);

    musp_i = mus_i * (1 - g);          % mu_s'
    mutp_i = mua_farrell + musp_i;     % mu_t'

    Ltot = sum(Li);
    if Ltot <= 0
        error('Sum of Li is zero/negative for k=%d.', k);
    end

    mutp_eff = Ltot / sum(Li ./ mutp_i);
    mu_sp_eff_all(idx) = mutp_eff - mua_farrell;   % effective mu_s'
end

R_Farrell = nan(numel(ks), numel(r_farrell));
for idx = 1:numel(ks)
    R_Farrell(idx,:) = dipole_reflectance_farrell(r_farrell, mua_farrell, mu_sp_eff_all(idx), n_farrell);
end

%% =========================================================
% (3) Schmitt layered model (REVERTED to earlier version)
%     Uses mu_eff = sqrt(3*mua*mus) with mus taken directly from sim_info (saved us)
% =========================================================
C_results = zeros(numel(ks),2);
R_Schmitt = nan(numel(ks), numel(rCenters));

for ii = 1:numel(ks)
    k = ks(ii);

    % Use MC radial curve already computed
    Rmc = R_MC(ii,:);

    % pathlengths
    Li = LayerPath{k}.Li(:);
    if numel(Li) < 4
        error('LayerPath{%d}.Li must have at least 4 elements for Schmitt.', k);
    end
    Li = Li(1:4);

    Ltop = sum(Li(1:3));
    if Ltop <= 0
        error('Sum of Li(1:3) is zero/negative for k=%d.', k);
    end

    % load sim_info optical props
    S = load(fullfile(files(k).folder, files(k).name), 'sim_info');
    if ~isfield(S,'sim_info')
        error('sim_info not found inside %s', fullfile(files(k).folder, files(k).name));
    end

    mua = double([S.sim_info.ua1, S.sim_info.ua2, S.sim_info.ua3, S.sim_info.ua4].');
    mus = double([S.sim_info.us1, S.sim_info.us2, S.sim_info.us3, S.sim_info.us4].');  % saved "us"

    % pathlength-weighted averages (top 3 collapsed)
    mua_top = sum(mua(1:3) .* Li(1:3)) / Ltop;
    mus_top = sum(mus(1:3) .* Li(1:3)) / Ltop;

    % bulk = layer 4
    mua_bulk = mua(4);
    mus_bulk = mus(4);

    % reverted mu_eff definition
    mu_eff1 = sqrt(3 * mua_top  * mus_top);
    mu_eff2 = sqrt(3 * mua_bulk * mus_bulk);

    % fit region
    valid = (rCenters > valid_r_min) & (rCenters < valid_r_max) & isfinite(Rmc) & (Rmc > 0);
    rfit  = rCenters(valid);
    yfit  = rfit .* Rmc(valid);     % r·R(r)

    schmittFun = @(C,rv) C(1)*exp(-mu_eff1*rv) + C(2)*exp(-mu_eff2*rv);
    C0 = [max(yfit), max(yfit)/5];

    if exist('lsqcurvefit','file') == 2
        opts = optimoptions('lsqcurvefit','Display','off');
        C = lsqcurvefit(schmittFun, C0, rfit, yfit, [], [], opts);
    else
        obj = @(C) sum((schmittFun(C,rfit) - yfit).^2);
        C = fminsearch(obj, C0);
    end

    C_results(ii,:) = C(:).';

    % prediction
    num = C(1)*exp(-mu_eff1*rCenters) + C(2)*exp(-mu_eff2*rCenters);
    Rmodel = num ./ rCenters;
    Rmodel(Rmodel <= 0) = NaN;  % for log plot safety

    R_Schmitt(ii,:) = Rmodel;
end

disp('Schmitt model coefficients [C1  C2] (reverted Schmitt):');
disp(array2table(C_results, 'VariableNames',{'C1','C2'}, 'RowNames',labels_sch));

%% =========================================================
% FIGURE A: All curves overlaid (single axes)
% =========================================================
figure('Color','w'); hold on;

plot(rCenters, R_MC(1,:), 'k-', 'LineWidth', 1.0, 'DisplayName','MC(model): No dye');
plot(rCenters, R_MC(2,:), 'r-', 'LineWidth', 1.0, 'DisplayName','MC(model): Dye 1');
plot(rCenters, R_MC(3,:), 'b-', 'LineWidth', 1.0, 'DisplayName','MC(model): Dye 2');

plot(r_farrell, R_Farrell(1,:), 'k--', 'LineWidth', 1.0, 'DisplayName','Farrell: No dye');
plot(r_farrell, R_Farrell(2,:), 'r--', 'LineWidth', 1.0, 'DisplayName','Farrell: Dye 1');
plot(r_farrell, R_Farrell(3,:), 'b--', 'LineWidth', 1.0, 'DisplayName','Farrell: Dye 2');

plot(rCenters, R_Schmitt(1,:), 'k-.', 'LineWidth', 1.0, 'DisplayName','Schmitt: No dye');
plot(rCenters, R_Schmitt(2,:), 'r-.', 'LineWidth', 1.0, 'DisplayName','Schmitt: Dye 1');
plot(rCenters, R_Schmitt(3,:), 'b-.', 'LineWidth', 1.0, 'DisplayName','Schmitt: Dye 2');

set(gca,'YScale','log');
grid on;
xlabel('Distance r (cm)');
ylabel('R (per area, per launch)');
title('Diffuse reflectance: MC(model-only) vs Farrell vs Schmitt');
legend('show','Location','northeast');
xlim(xlim_plot);
ylim(ylim_plot);
set(gca,'FontName','Arial','FontSize',10);

%% =========================================================
% FIGURE B: Subplots (MC / Farrell / Schmitt)
% =========================================================
figure('Color','w','Position',[100 100 1200 320]);

subplot(1,3,1); hold on;
plot(rCenters, R_MC(1,:), 'k-', 'LineWidth', 1.0, 'DisplayName','No dye');
plot(rCenters, R_MC(2,:), 'r-', 'LineWidth', 1.0, 'DisplayName','1 layer dye');
plot(rCenters, R_MC(3,:), 'b-', 'LineWidth', 1.0, 'DisplayName','2 layer dye');
set(gca,'YScale','log'); grid on; axis square;
xlabel('r (cm)'); ylabel('R (a.u)'); title('(a) Monte Carlo');
xlim(xlim_plot); ylim(ylim_plot); legend('Location','southwest');

subplot(1,3,2); hold on;
plot(r_farrell, R_Farrell(1,:), 'k-', 'LineWidth', 1.0, 'DisplayName','No dye');
plot(r_farrell, R_Farrell(2,:), 'r-', 'LineWidth', 1.0, 'DisplayName','1 layer dye');
plot(r_farrell, R_Farrell(3,:), 'b-', 'LineWidth', 1.0, 'DisplayName','2 layer dye');
set(gca,'YScale','log'); grid on; axis square;
xlabel('r (cm)'); ylabel('R (a.u)'); title('(b) Farrell dipole');
xlim(xlim_plot); ylim(ylim_plot); legend('Location','southwest');

subplot(1,3,3); hold on;
plot(rCenters, R_Schmitt(1,:), 'k-', 'LineWidth', 1.0, 'DisplayName','No dye');
plot(rCenters, R_Schmitt(2,:), 'r-', 'LineWidth', 1.0, 'DisplayName','1 layer dye');
plot(rCenters, R_Schmitt(3,:), 'b-', 'LineWidth', 1.0, 'DisplayName','2 layer dye');
set(gca,'YScale','log'); grid on; axis square;
xlabel('r (cm)'); ylabel('R (a.u)'); title('(c) Schmitt layered');
xlim(xlim_plot); ylim(ylim_plot); legend('Location','southwest');

set(findall(gcf,'-property','FontName'),'FontName','Arial');
set(findall(gcf,'-property','FontSize'),'FontSize',10);

%% =========================================================
% FIGURE C: 2x2 comparison (MC / Exp / Farrell / Schmitt)
% Layout:
%   (1,1) MC        (1,2) EXP
%   (2,1) Farrell   (2,2) Schmitt
%% =========================================================

% ---- Experimental data (baseline voltages) ----
distance_cm = [2 3 4.5];
V_noDye_mV  = [0.463  0.111  0.022];
V_dye_mV    = [0.6018 0.173  0.0413];

figure('Color','w','Position',[100 100 900 700]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% ---------------- (1,1) MC ----------------
nexttile; hold on; box on; grid on;
plot(rCenters, R_MC(1,:), 'k-', 'LineWidth', 1, 'DisplayName','No dye');
plot(rCenters, R_MC(2,:), 'r-', 'LineWidth', 1, 'DisplayName','1 layer dye');
plot(rCenters, R_MC(3,:), 'b-', 'LineWidth', 1, 'DisplayName','2 layer dye');
set(gca,'YScale','log');
xlabel('r (cm)'); ylabel('R (a.u.)');
title('(a) MC model');
xlim([0 6]);
ylim([1e-5 10]);
legend('Location','southwest');
axis square;

% ---------------- (1,2) EXP ----------------
nexttile; hold on; box on; grid on;
plot(distance_cm, V_noDye_mV, 'ko-', ...
    'LineWidth', 1, 'MarkerSize', 7, 'DisplayName','No dye');
plot(distance_cm, V_dye_mV, 'ro-', ...
    'LineWidth', 1, 'MarkerSize', 7, 'DisplayName','1 layer dye');
xlabel('r (cm)');
ylabel('R (mV)');
title('(b) Experiment (chicken)');
xlim([0 6]);
legend('Location','northeast');
axis square;

% ---------------- (2,1) Farrell ----------------
nexttile; hold on; box on; grid on;
plot(r_farrell, R_Farrell(1,:), 'k-', 'LineWidth', 1, 'DisplayName','No dye');
plot(r_farrell, R_Farrell(2,:), 'r-', 'LineWidth', 1, 'DisplayName','1 layer dye');
plot(r_farrell, R_Farrell(3,:), 'b-', 'LineWidth', 1, 'DisplayName','2 layer dye');
set(gca,'YScale','log');
xlabel('r (cm)'); ylabel('R (a.u.)');
title('(c) Farrell dipole');
xlim([0 6]);
ylim([1e-5 10]);
legend('Location','southwest');
axis square;

% ---------------- (2,2) Schmitt ----------------
nexttile; hold on; box on; grid on;
plot(rCenters, R_Schmitt(1,:), 'k-', 'LineWidth', 1, 'DisplayName','No dye');
plot(rCenters, R_Schmitt(2,:), 'r-', 'LineWidth', 1, 'DisplayName','1 layer dye');
plot(rCenters, R_Schmitt(3,:), 'b-', 'LineWidth', 1, 'DisplayName','2 layer dye');
set(gca,'YScale','log');
xlabel('r (cm)'); ylabel('R (a.u.)');
title('(d) Schmitt layered');
xlim([0 6]);
ylim([1e-5 10]);
legend('Location','southwest');
axis square;

% ---- Global formatting ----
set(findall(gcf,'-property','FontName'),'FontName','Arial');
set(findall(gcf,'-property','FontSize'),'FontSize',10);




%% =========================================================
% Local function: Farrell dipole
% =========================================================
function R = dipole_reflectance_farrell(r, mua, musp, n)
r = r(:)'; % row vector

mutp  = mua + musp;                % mu_t' [cm^-1]
D     = 1 ./ (3 * mutp);           % diffusion coefficient [cm]
mueff = sqrt(mua ./ D);            % [cm^-1]

Reff = -1.440./n.^2 + 0.710./n + 0.668 + 0.0636.*n;
A    = (1 + Reff) ./ (1 - Reff);

z0 = 1 ./ mutp;                    % [cm]
zb = 2 .* A .* D;                  % [cm]

r1 = sqrt(r.^2 + z0.^2);
r2 = sqrt(r.^2 + (z0 + 2*zb).^2);

R = (1/(4*pi)) .* ( ...
      z0 .* (mueff + 1./r1) .* exp(-mueff .* r1) ./ (r1.^2) + ...
     (z0 + 2*zb) .* (mueff + 1./r2) .* exp(-mueff .* r2) ./ (r2.^2) );
end

%% =========================================================
% Local function: radial average on LINEAR normalized map
% =========================================================
function [rCenters, Ravg] = radialAverage_linMap(Rmap_lin, xEdges, yEdges)

xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;

[X,Y] = meshgrid(xCenters, yCenters);
Rrad  = sqrt(X.^2 + Y.^2);

Rmap = Rmap_lin.';  % match X,Y (Ny×Nx)

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
end
