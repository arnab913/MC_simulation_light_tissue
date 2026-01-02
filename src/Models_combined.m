%% =========================================================
%  ALL MODELS IN ONE SCRIPT:
%   (1) MC radial reflectance
%   (2) Farrell diffusion dipole
%   (3) Schmitt layered model
%  PLUS: a final figure with 3 subplots (MC / Farrell / Schmitt)
% =========================================================

% ---- housekeeping (comment out if you don't want clearing) ----
% clearvars -except Diff_refl xEdges yEdges LayerPath files
% close all; clc;

%% -----------------------------
% Shared settings
% -----------------------------
ks = [1 6 11];

% Radial binning from MC grid (cm)
xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
[X,Y]    = meshgrid(xCenters, yCenters);
Rrad     = sqrt(X.^2 + Y.^2);

dr       = mean(diff(xCenters));
rBins    = 0:dr:max(xCenters);
rCenters = (rBins(1:end-1) + rBins(2:end)) / 2;

% Farrell dipole settings
mua_farrell = 0.05;                   % [cm^-1]
n_farrell   = 1.4;                    % refractive index
r_farrell   = linspace(0.01, 6, 400); % [cm]

% Farrell layer reduced scattering [cm^-1] (EDIT if needed)
musp_layers = struct();
musp_layers.k1  = [30 30 30 30];
musp_layers.k6  = [5 10 20 30];
musp_layers.k11 = [2 5 10 30];

% Schmitt settings
labels_sch  = {'No dye','1 layer dye','2 layer dye'};
valid_r_min = 0.3;
valid_r_max = 4.0;

% =========================================================
% (1) MC radial reflectance (compute & store)
% =========================================================
R_MC = nan(numel(ks), numel(rCenters));  % rows correspond to ks

for a = 1:numel(ks)
    kk = ks(a);

    if numel(Diff_refl) < kk || isempty(Diff_refl{kk}) || ~isfield(Diff_refl{kk},'map')
        error('Diff_refl{%d}.map not found.', kk);
    end

    RRlog = Diff_refl{kk}.map;  % log10(R)
    RRlin = 10.^RRlog;          % linear R

    R_mean = nan(size(rCenters));
    for i = 1:numel(rCenters)
        mask = (Rrad >= rBins(i)) & (Rrad < rBins(i+1));
        vals = RRlin(mask);
        if ~isempty(vals)
            R_mean(i) = mean(vals,'omitnan');
        end
    end

    R_MC(a,:) = R_mean;
end

% =========================================================
% (2) Farrell diffusion dipole (compute & store)
% =========================================================
mu_sp_eff_all = zeros(size(ks));

for idx = 1:numel(ks)
    k = ks(idx);

    if numel(LayerPath) < k || isempty(LayerPath{k}) || ~isfield(LayerPath{k},'Li')
        error('LayerPath{%d}.Li not found. Make sure LayerPath is in workspace.', k);
    end

    Li = LayerPath{k}.Li;
    if numel(Li) < 4
        error('LayerPath{%d}.Li must have at least 4 elements.', k);
    end

    Li = Li(1:4);
    Li = Li(:); % column
    if any(~isfinite(Li)) || any(Li < 0)
        error('LayerPath{%d}.Li contains invalid values.', k);
    end

    switch k
        case 1
            musp_i = musp_layers.k1;
        case 6
            musp_i = musp_layers.k6;
        case 11
            musp_i = musp_layers.k11;
        otherwise
            error('Unexpected k=%d', k);
    end
    musp_i = musp_i(:);

    Ltot = sum(Li);
    if Ltot == 0
        error('Sum of Li is zero for k=%d.', k);
    end

    mutp_i   = mua_farrell + musp_i;           % mu_t'
    mutp_eff = Ltot / sum(Li ./ mutp_i);       % harmonic mean along path
    mu_sp_eff_all(idx) = mutp_eff - mua_farrell;
end

R_Farrell = nan(numel(ks), numel(r_farrell));
for idx = 1:numel(ks)
    R_Farrell(idx,:) = dipole_reflectance_farrell(r_farrell, mua_farrell, mu_sp_eff_all(idx), n_farrell);
end

% =========================================================
% (3) Schmitt layered model (fit & store)
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
    if Ltop == 0
        error('Sum of Li(1:3) is zero for k=%d.', k);
    end

    % load sim_info optical props
    if ~exist('files','var') || numel(files) < k
        error('files struct not found or files(%d) missing. Needed to load sim_info for Schmitt.', k);
    end
    S = load(fullfile(files(k).folder, files(k).name), 'sim_info');
    if ~isfield(S,'sim_info')
        error('sim_info not found inside %s', fullfile(files(k).folder, files(k).name));
    end

    mua = double([S.sim_info.ua1, S.sim_info.ua2, S.sim_info.ua3, S.sim_info.ua4].');
    mus = double([S.sim_info.us1, S.sim_info.us2, S.sim_info.us3, S.sim_info.us4].');

    % pathlength-weighted averages for top (layers 1–3)
    mua_top = sum(mua(1:3) .* Li(1:3)) / Ltop;
    mus_top = sum(mus(1:3) .* Li(1:3)) / Ltop;

    % bulk = layer 4
    mua_bulk = mua(4);
    mus_bulk = mus(4);

    mu_eff1 = sqrt(3 * mua_top  * mus_top);
    mu_eff2 = sqrt(3 * mua_bulk * mus_bulk);

    valid = (rCenters > valid_r_min) & (rCenters < valid_r_max) & isfinite(Rmc) & (Rmc > 0);
    rfit  = rCenters(valid);
    yfit  = rfit .* Rmc(valid);

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

    Rmodel = (C(1)*exp(-mu_eff1*rCenters) + C(2)*exp(-mu_eff2*rCenters)) ./ rCenters;
    R_Schmitt(ii,:) = Rmodel;
end

disp('Schmitt model coefficients [C1  C2]:');
disp(array2table(C_results, 'VariableNames',{'C1','C2'}, 'RowNames',labels_sch));

%% =========================================================
% FIGURE A: All curves overlaid (single axes)
% =========================================================
figure('Color','w'); hold on;

for i = 1:numel(ks)
    plot(rCenters, R_MC(i,:), 'LineWidth', 1.0, ...
        'DisplayName', sprintf('MC: k=%d', ks(i)));
end
plot(r_farrell, R_Farrell(1,:), 'k-', 'LineWidth', 1.5, 'DisplayName','Farrell: No dye');
plot(r_farrell, R_Farrell(2,:), 'r-', 'LineWidth', 1.5, 'DisplayName','Farrell: Dye 1');
plot(r_farrell, R_Farrell(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName','Farrell: Dye 2');

for i = 1:numel(ks)
    plot(rCenters, R_Schmitt(i,:), '-', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Schmitt: %s', labels_sch{i}));
end

set(gca,'YScale','log');
grid on;
xlabel('Distance r (cm)');
ylabel('R (a.u.)');
title('Diffuse reflectance: MC vs Farrell vs Schmitt (overlay)');
legend('show','Location','northeast');
xlim([0 6]);
set(gca,'FontName','Arial','FontSize',10);

%% =========================================================
% FIGURE B: Subplots (MC / Farrell / Schmitt)
% =========================================================
figure('Color','w');

% --- (B1) MC subplot ---
subplot(1,3,1); hold on;
for i = 1:numel(ks)
    plot(rCenters, R_MC(i,:), 'LineWidth', 1, ...
        'DisplayName', sprintf('k=%d', ks(i)));
end
set(gca,'YScale','log'); grid on;
xlabel('r (cm)'); ylabel('R (a.u.)');
axis square;
title('(a) MC');
xlim([0 6]);
set(gca,'FontName','Arial','FontSize',10);
legend('show','Location','southwest');

% --- (B2) Farrell subplot ---
subplot(1,3,2); hold on;
plot(r_farrell, R_Farrell(1,:), 'k-', 'LineWidth', 1, 'DisplayName','No dye');
plot(r_farrell, R_Farrell(2,:), 'r-', 'LineWidth', 1, 'DisplayName','Dye 1');
plot(r_farrell, R_Farrell(3,:), 'b-', 'LineWidth', 1, 'DisplayName','Dye 2');
set(gca,'YScale','log'); grid on;
xlabel('r (cm)'); ylabel('R (a.u.)');
title('(b) Farrell dipole');
axis square;
xlim([0 6]);
set(gca,'FontName','Arial','FontSize',10);
legend('show','Location','southwest');

% --- (B3) Schmitt subplot ---
subplot(1,3,3); hold on;
for i = 1:numel(ks)
    plot(rCenters, R_Schmitt(i,:), 'LineWidth', 1, ...
        'DisplayName', labels_sch{i});
end
set(gca,'YScale','log'); grid on;
xlabel('r (cm)'); ylabel('R (a.u.)');
title('(c) Schmitt layered');
axis square;
xlim([0 6]);
set(gca,'FontName','Arial','FontSize',10);
legend('show','Location','southwest');

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
