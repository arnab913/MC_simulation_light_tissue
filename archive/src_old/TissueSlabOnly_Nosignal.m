% --- Sweep parameters ---
ua_ink = [1,7];                   % for g= 0.9, us is 10 times us' 
results_all = zeros(1,7); % Reflectance vs distance (for 1-7 cm)
all_det_r_cm = 1:0.5:6;
for idx = 1:length(ua_ink)
    % --- Model and geometry setup ---
    model = MCmatlab.model;
    model.G.nx = 101;
    model.G.ny = 101;
    model.G.nz = 150;
    model.G.Lx = 15;    % [cm]
    model.G.Ly = 15;
    model.G.Lz = 5;
   %model.MC.nPhotons = 100000;
    model.G.mediaPropertiesFunc = @mediaPropertiesFunc;
    model.G.geomFunc = @geometryDefinition;
    model.G.mediaPropParams = {ua_ink(idx)}; % Pass scattering coeff as parameter

    % --- MC simulation parameters ---
  %  model.MC.simulationTimeRequested = 2;   % [min]
    model.MC.nPhotonsRequested = 5e7;
    model.MC.matchedInterfaces = true;
    model.MC.boundaryType = 1;
    model.MC.wavelength = 760;                % [nm]
    %--- LED source (simple center, top-hat spatial, Lambertian angular) ---
    model.MC.lightSource.sourceType = 0;
    model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 1;
    model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 1;
    model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = 0.1; % [cm]
    model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = 0.1; % [cm]
    model.MC.lightSource.angularIntensityDistribution.XDistr = 2; % Lambertian
    model.MC.lightSource.angularIntensityDistribution.YDistr = 2; % Lambertian
    model.MC.lightSource.theta = 0;
    model.MC.lightSource.phi = 0;
    model.MC.lightSource.psi = 0;
    % Place source at center surface:
    model.MC.lightSource.xFocus =  0; % Center
    model.MC.lightSource.yFocus =  0; % Center
    model.MC.lightSource.zFocus =  0; % Just at surface

    % --- Run simulation ---
    model = plot(model,'G');
    model = runMonteCarlo(model);
    model = plot(model,'MC');


    % --- Extract reflectance at various distances ---
    nx = model.G.nx; ny = model.G.ny; Lx = model.G.Lx; Ly = model.G.Ly;
    x = linspace(-Lx/2, Lx/2, nx);
    y = linspace(-Ly/2, Ly/2, ny);
    [X, Y] = meshgrid(x, y);  % 2D grid for surface

    source_x = 0;
    source_y = 0;
    dist_from_source = sqrt((X - source_x).^2 + (Y - source_y).^2);
    Rmap = model.MC.normalizedIrradiance_zneg;
    det_radius = 0.3; % [cm]
    reflectance = zeros(size(all_det_r_cm));
 
  for j = 1:length(all_det_r_cm)
       % Detector center position
       det_x = all_det_r_cm(j);
       det_y = 0; % same y as source
    
       % Create a circular mask for detector area
       mask = ((X - det_x).^2 + (Y - det_y).^2) <= det_radius^2;
    
       % Sum (or mean) reflectance within detector area
       reflectance(j) = sum(Rmap(mask));
  end

    % --- Save results for this run ---
    % fname = sprintf('TestTissue_us %d Baseline.mat',ua_ink(idx));
    % save(fname, 'reflectance', 'all_det_r_cm');,
    fname_model = sprintf('ND_Baseline_%g.mat', ua_ink(idx));
    save(fname_model, 'model');  % saves model.G, model.MC, etc.

    % (Optional) Display progress
    fprintf('Finished us = %d cm^-1\n',ua_ink(idx));
end

% Main simulation script here
... (all your setup code) ...

% === Supporting Functions (must be at end of this file) ===
function M = geometryDefinition(X,Y,Z,parameters)
    x_cyl_min = -7.5; % Cylinder runs from -Lx/2 to Lx/2 (full extent)
    x_cyl_max = 7.5;
    y_cyl = 0;        % Center in y at 0
    z_cyl = 1;      % Center in z at 1.0 cm
    r_cyl = 0.3;      % Radius 0.3 cm (0.6 cm diameter)

    % Tissue everywhere by default
    M = ones(size(X));

    % Horizontal cylinder: axis along x, so
    % Condition: 
    % - Within x extent, AND
    % - Cross-sectional circle in Y-Z: (Y-y_cyl)^2 + (Z-z_cyl)^2 <= r^2
    % mask = (X >= x_cyl_min) & (X <= x_cyl_max) & ...
    %        ((Y - y_cyl).^2 + (Z - z_cyl).^2 <= r_cyl^2);
    % 
    % M(mask) = 2;   % Media 2 is the ink
end


function mediaProperties = mediaPropertiesFunc(param)
    mediaProperties = MCmatlab.mediumProperties;
    j=1;
    mediaProperties(j).name = 'Tissue';
    mediaProperties(j).mua = 0.05;
    mediaProperties(j).mus = 30;
    mediaProperties(j).g   = 0.9;

    %Biomed Opt Express. 2014 Jun 4;5(7):2037–2053. doi: 10.1364/BOE.5.002037..
    j=2;
    mediaProperties(j).name = 'Ink layer';
    mediaProperties(j).mua = param{1};   % mimicking 10% ink
    mediaProperties(j).mus = 0.8;   % also check Zotero/MC simulation part
    mediaProperties(j).g   = 0.8;
end


%% plot all

ua_ink = [1,7];
color_list = [
   0.85 0.33 0.10;   % orange
   0.13 0.65 0.14;   % green
   0.14 0.44 0.82;   % blue
   0.75 0.15 0.55;   % purple
   0.55 0.52 0.09    % olive
]; % 5 unique colors (rows are RGB)

files = dir('*.mat');
figure;
legendLabels = cell(1, numel(files));
for k = 1:numel(files)
    data = load(files(k).name);
    splitnum = regexp(files(k).name, '\d+', 'match');
    if ~isempty(splitnum)
        us = str2double(splitnum{1}); % Assumes first number is us
        color_idx = find(ua_ink == us, 1);
        if isempty(color_idx)
            plot_color = [0.5 0.5 0.5]; % if not found, use gray
        else
            plot_color = color_list(color_idx, :);
        end
    else
        plot_color = [0.5 0.5 0.5]; % fallback gray
    end
    % Solid for ink, dashed for baseline
    if contains(files(k).name, 'Baseline', 'IgnoreCase', true)
        style = 's--';
    else
        style = 'o-';
    end
    semilogy(data.all_det_r_cm, data.reflectance, style, 'Color', plot_color, 'LineWidth', 1.8);
    hold on;
    legendLabels{k} = files(k).name;
end
hold off;
xlabel('Source–detector distance (cm)');
ylabel('Intensity (a.u.)');
title('Baseline and Ink signal Only Tissue');
legend(legendLabels, 'Location', 'best');
grid on;
setJournalFigure(gcf, 'Interpreter', 'tex');