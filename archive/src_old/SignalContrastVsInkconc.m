% load all files
ND_base=load("ND_Baseline_1.mat");
ND_1cm = load("ND_ink_1.mat");
ND_7cm = load("ND_ink_7.mat");
ND_30cm = load("ND_ink_30.mat");
ND_100cm =load("ND_ink_100.mat");

D1_base= load("Dye_Baseline_1.mat");
D1_1cm = load("Dye_ink_1.mat");
D1_7cm = load("Dye_ink_7.mat");
D1_30cm = load("Dye_ink_30.mat");
D1_100cm = load("Dye_ink_100.mat");

%% Common geometry info (take from any loaded model)
model_ref = ND_base.model;
nx = model_ref.G.nx; ny = model_ref.G.ny;
Lx = model_ref.G.Lx; Ly = model_ref.G.Ly;

x = linspace(-Lx/2, Lx/2, nx);
y = linspace(-Ly/2, Ly/2, ny);
[X, Y] = meshgrid(x, y);

%% Detector geometry: 6 cm from source along +x, radius 0.3 cm
source_x = 0; source_y = 0;
det_r_cm = 6;             % detector center distance from source
det_radius = 0.3;         % [cm]
det_x = det_r_cm;         % center position
det_y = 0;

det_mask = ((X - det_x).^2 + (Y - det_y).^2) <= det_radius^2;

%% Helper function: integrate reflectance in detector area
get_Rdet = @(M) sum(M.MC.normalizedIrradiance_zneg(det_mask));

%% Ink absorption values [cm^-1] including baseline as 0
%ua_ink = [0 1 7 30 100];
ua_ink = [0 0.05 0.25 1 5];    %in percentage

%% Compute detector signal for ND (no dye layers)
R_ND_base  = get_Rdet(ND_base.model);
R_ND_1     = get_Rdet(ND_1cm.model);
R_ND_7     = get_Rdet(ND_7cm.model);
R_ND_30    = get_Rdet(ND_30cm.model);
R_ND_100   = get_Rdet(ND_100cm.model);
R_ND_all   = [R_ND_base, R_ND_1, R_ND_7, R_ND_30, R_ND_100];

%% Compute detector signal for Dye Layer 1
R_D1_base  = get_Rdet(D1_base.model);
R_D1_1     = get_Rdet(D1_1cm.model);
R_D1_7     = get_Rdet(D1_7cm.model);
R_D1_30    = get_Rdet(D1_30cm.model);
R_D1_100   = get_Rdet(D1_100cm.model);
R_D1_all   = [R_D1_base, R_D1_1, R_D1_7, R_D1_30, R_D1_100];

%% Compute detector signal for Dye Layer 2
R_D2_base  = get_Rdet(D2_base.model);
R_D2_1     = get_Rdet(D2_1cm.model);
R_D2_7     = get_Rdet(D2_7cm.model);
R_D2_30    = get_Rdet(D2_30cm.model);
R_D2_100   = get_Rdet(D2_100cm.model);
R_D2_all   = [R_D2_base, R_D2_1, R_D2_7, R_D2_30, R_D2_100];

%% Compute contrast = (with ink - corresponding baseline)
% Baseline is the 0-th ink entry (index 1)
C_ND = R_ND_all - R_ND_base;  % ND case (may just be 0 for baseline)
C_D1 = R_D1_all - R_D1_base;  % Dye layer 1
C_D2 = R_D2_all - R_D2_base;  % Dye layer 2

% We usually don't plot the baseline point (0 contrast at ua_ink=0),
% but you can if you want. Here we plot only ink concentrations >0:
ua_ink_wo0 = ua_ink(2:end);
C_ND_wo0   = C_ND(2:end);
C_D1_wo0   = C_D1(2:end);
C_D2_wo0   = C_D2(2:end);

%% Plot: contrast vs ink absorption at 6 cm, 0.3 cm radius detector
figure;
semilogx(ua_ink_wo0, C_ND_wo0, 'o-',  'LineWidth',1.5); hold on;
semilogx(ua_ink_wo0, C_D1_wo0, 's--', 'LineWidth',1.5);
semilogx(ua_ink_wo0, C_D2_wo0, 'd-.', 'LineWidth',1.5);
hold off;

xlabel('\mu_a^{ink} (cm^{-1})');
ylabel('Contrast at 6 cm (ink - baseline) [a.u.]');
title('Signal Contrast vs Ink Absorption (6 cm, r = 0.3 cm)');
legend({'ND (no dye)','Dye Layer 1','Dye Layer 2'}, 'Location','best');
grid on;
set(gca,'XScale','log');  % log scale over 1–100 cm^-1
