% load all files
ND_base=load("ND_base_Model.mat");
ND_1cm = load("MCmodel_ua_1.mat");
ND_7cm = load("MCmodel_ua_7.mat");
ND_30cm = load("MCmodel_ua_30.mat");
ND_100cm =load("MCmodel_ua_100.mat");

D1_base= load("DyeL1_base_Model.mat");
D1_1cm = load("DyeL1_MCmodel_ink_1.mat");
D1_7cm = load("DyeL1_MCmodel_ink_7.mat");
D1_30cm = load("DyeL1_MCmodel_ink_30.mat");
D1_100cm = load("DyeL1_MCmodel_ink_100.mat");

D2_base= load("DyeL2_Base_MCmodel.mat");
D2_1cm = load("DyeL2_MCmodel_ink_1.mat");
D2_7cm = load("DyeL2_MCmodel_ink_7.mat");
D2_30cm = load("DyeL2_MCmodel_ink_30.mat");
D2_100cm = load("DyeL2_MCmodel_ink_100.mat");



% Extract middle-slice (y = 51) for each model
ND = rot90(log(squeeze(ND_base.model.MC.normalizedFluenceRate(:,51,:))));
D1 = rot90(log(squeeze(D1_base.model.MC.normalizedFluenceRate(:,51,:))));
D2 = rot90(log(squeeze(D2_base.model.MC.normalizedFluenceRate(:,51,:))));

% Find global min/max for same color scale
allVals = [ND(:); D1(:); D2(:)];
cmin = min(allVals);
cmax = max(allVals);


figure;

subplot(1,3,1);
contourf(ND, 50);
title('ND (No Dye)');
caxis([cmin cmax]);
ylim([100 150])
axis square
colorbar;

subplot(1,3,2);
contourf(D1, 50);
title('D1 (Dye Layer 1)');
caxis([cmin cmax]);
ylim([100 150])
axis square
colorbar;

subplot(1,3,3);
contourf(D2, 50);
title('D2 (Dye Layer 2)');
caxis([cmin cmax]);
ylim([100 150])
axis square
colorbar;
tightfig; 
sgtitle('XZ Fluence Maps (same colorscale)');


%%
% Extract the three matrices first
A = rot90(log(squeeze(ND_base.model.MC.normalizedFluenceRate(:,:,1))));
B = rot90(log(squeeze(D1_base.model.MC.normalizedFluenceRate(:,:,1))));
C = rot90(log(squeeze(D2_base.model.MC.normalizedFluenceRate(:,:,1))));

% Find universal color scale
cmin = min([A(:); B(:); C(:)]);
cmax = max([A(:); B(:); C(:)]);

figure;
sgtitle('Diffuse Reflectance at Top')

% ---- Subplot 1 ----
subplot(1,3,1)
contourf(A)
axis square
%ylim([100 150])
caxis([cmin cmax])
title('ND')
colorbar

% ---- Subplot 2 ----
subplot(1,3,2)
contourf(B)
axis square
%ylim([100 150])
caxis([cmin cmax])
title('D1')
colorbar

% ---- Subplot 3 ----
subplot(1,3,3)
contourf(A)
axis square
%ylim([100 150])
caxis([cmin cmax])
title('D2')
colorbar

%% 

% Find universal color scale
cmin = min([B(:); C(:)]);
cmax = max([B(:); C(:)]);
figure()
surf((rot90(squeeze(D2_base.model.MC.normalizedFluenceRate(51-30:51+30,51-30:51+30,1)-squeeze(ND_base.model.MC.normalizedFluenceRate(51-30:51+30,51-30:51+30,1))))));axis square
%set(gca, 'ColorScale', 'log');

figure;
sgtitle('Diffuse Reflectance Increase after dye');

% ---- Subplot 1 ----
subplot(1,2,1)
diff_BA =B-A;
surf(diff_BA);
axis square
%ylim([100 150])
caxis([cmin cmax])
title('Baseline difference (Dye Layer 1-No Dye)');
colorbar

% ---- Subplot 2 ----
subplot(1,2,2)
diff_CA =C-A;
surf(diff_CA);
axis square
%ylim([100 150])
caxis([cmin cmax])
title('Baseline difference (Dye Layer 2-No Dye)');
colorbar

%% Photon increase at ink layers

% ---- Common settings ----
y_idx   = 51;                      % middle slice in y
nLevels = 50;                      % contour levels
inkVals = [1 7 30 100];            % ink concentrations (cm^-1)

% Pack models into cell arrays for looping
ND_models = {ND_1cm, ND_7cm, ND_30cm, ND_100cm};
D1_models = {D1_1cm, D1_7cm, D1_30cm, D1_100cm};
D2_models = {D2_1cm, D2_7cm, D2_30cm, D2_100cm};

for ii = 1:numel(inkVals)

    % ---- Extract XZ log-fluence slices for this ink value ----
    ND_slice = rot90(log(squeeze(ND_models{ii}.model.MC.normalizedFluenceRate(:,y_idx,:))));
    D1_slice = rot90(log(squeeze(D1_models{ii}.model.MC.normalizedFluenceRate(:,y_idx,:))));
    D2_slice = rot90(log(squeeze(D2_models{ii}.model.MC.normalizedFluenceRate(:,y_idx,:))));

    % Compute global color scale for the three maps of this ink conc.
    allVals = [ND_slice(:); D1_slice(:); D2_slice(:)];
    cmin = min(allVals);
    cmax = max(allVals);

    % ---- Plot 1x3 subplots ----
    figure;

    subplot(1,3,1);
    contourf(ND_slice, nLevels);
    title(sprintf('Depth 1 cm, Ink %g cm^{-1} (No dye)', inkVals(ii)));
    caxis([cmin cmax]); axis square; colorbar;

    subplot(1,3,2);
    contourf(D1_slice, nLevels);
    title(sprintf('Depth 1 cm, Ink %g cm^{-1} (Dye Layer 1)', inkVals(ii)));
    caxis([cmin cmax]); axis square; colorbar;

    subplot(1,3,3);
    contourf(D2_slice, nLevels);
    title(sprintf('Depth 1 cm, Ink %g cm^{-1} (Dye Layer 2)', inkVals(ii)));
    caxis([cmin cmax]); axis square; colorbar;

    sgtitle(sprintf('XZ Fluence Maps, Ink %g cm^{-1} (same colorscale)', inkVals(ii)));
end

%% calculate absorption at ink position 1cm depth
z_idx = 30;        % xy plane at this depth
nLevels = 50;

models_ND = {ND_1cm, ND_7cm, ND_30cm, ND_100cm};
models_D1 = {D1_1cm, D1_7cm, D1_30cm, D1_100cm};
models_D2 = {D2_1cm, D2_7cm, D2_30cm, D2_100cm};
inkVals   = [1 7 30 100];

for ii = 1:numel(inkVals)

    % --- Extract xy slices of normalizedAbsorption at z = z_idx ---
    A_ND = squeeze(models_ND{ii}.model.MC.normalizedAbsorption(:,:,z_idx));
    A_D1 = squeeze(models_D1{ii}.model.MC.normalizedAbsorption(:,:,z_idx));
    A_D2 = squeeze(models_D2{ii}.model.MC.normalizedAbsorption(:,:,z_idx));

    % Optional: log-scale for dynamic range
    A_ND_log = log(A_ND);
    A_D1_log = log(A_D1);
    A_D2_log = log(A_D2);

    % Common color scale per ink value
    allVals = [A_ND_log(:); A_D1_log(:); A_D2_log(:)];
    cmin = min(allVals);
    cmax = max(allVals);

    figure;

    subplot(1,3,1);
    contourf(A_ND_log, nLevels);
    title(sprintf('ND, ink %g cm^{-1}, z-idx %d', inkVals(ii), z_idx));
    caxis([cmin cmax]); axis square; colorbar;

    subplot(1,3,2);
    contourf(A_D1_log, nLevels);
    title(sprintf('Dye L1, ink %g cm^{-1}, z-idx %d', inkVals(ii), z_idx));
    caxis([cmin cmax]); axis square; colorbar;

    subplot(1,3,3);
    contourf(A_D2_log, nLevels);
    title(sprintf('Dye L2, ink %g cm^{-1}, z-idx %d', inkVals(ii), z_idx));
    caxis([cmin cmax]); axis square; colorbar;

    sgtitle(sprintf('XY normalizedAbsorption at ~1 cm depth (z = %d)', z_idx));
end
