function [x_in, y_in, z_in, x_ot, y_ot, z_ot, s, w, no_of_photons, M_raw] = ...
    common_tissue_simulation(ua1, us1, n1, ua2, us2, n2, ua3, us3, n3)

addpath(genpath('/home/kniknam/MCmatlab 4.4.9.0'))

%% Simulation parameters
nBatches       = 20;        % number of batches
nPhotonsReq    = 100000;      % photons per batch
nExamplePaths  = 5000;      % example paths per batch
cmd_size       = [10, 10, 5]; % simulation volume [cm]
zSurface       = 0;         % surface depth [cm]
wvlngth        = 760;       % wavelength [nm]
beam_X         = 0; beam_Y = 0; beam_phi = 0; beam_tht = 0;
dl             = 0.05;      % spatial resolution [cm]

%% Initialize MC model
model = MCmatlab.model;
model.G.nx = round(cmd_size(1)/dl);
model.G.ny = round(cmd_size(2)/dl);
model.G.nz = round(cmd_size(3)/dl);
model.G.Lx = cmd_size(1); model.G.Ly = cmd_size(2); model.G.Lz = cmd_size(3);

% Media properties
model.G.mediaPropertiesFunc = @mediaPropertiesFunc;
model.G.mediaPropParams = {ua1, us1, n1, ua2, us2, n2, ua3, us3, n3};
model.G.geomFunc = @geometryDefinition;

% Monte Carlo setup
model.MC.matchedInterfaces = false;
model.MC.boundaryType = 1;
model.MC.nPhotonsRequested = nPhotonsReq;
model.MC.nExamplePaths = nExamplePaths;
model.MC.wavelength = wvlngth;

% Pencil beam
model.MC.lightSource.sourceType = 0;
model.MC.silentMode = false;
model.MC.lightSource.xFocus = beam_X;
model.MC.lightSource.yFocus = beam_Y;
model.MC.lightSource.zFocus = zSurface;
model.MC.lightSource.theta  = beam_tht;
model.MC.lightSource.phi    = beam_phi;

scatter_stat = [];   % ALL batches combined

for batch = 1:nBatches
    disp(['Running batch ', num2str(batch), '/', num2str(nBatches)])

    batch_scatter_stat = [];   % RESET per batch

    while size(batch_scatter_stat,1) < nExamplePaths

        t_model = runMonteCarlo(model);
        paths = t_model.MC.examplePaths;
        spratrs = find(isnan(paths(1,:)));

        t_scatter_stat = nan(length(spratrs)-1, 8);
        cnt = 0;

        for i = 1:length(spratrs)-1
            start_col = spratrs(i)+1;
            end_col   = spratrs(i+1)-1;
            if start_col > end_col
                continue
            end

            Photon_Path = paths(:, start_col:end_col);
            cnt = cnt + 1;

            step_vecs = Photon_Path(1:3,2:end) - Photon_Path(1:3,1:end-1);
            step_len  = sqrt(sum(step_vecs.^2,1));

            t_scatter_stat(cnt,:) = [ ...
                Photon_Path(1,1), Photon_Path(2,1), Photon_Path(3,1), ...
                Photon_Path(1,end), Photon_Path(2,end), Photon_Path(3,end), ...
                sum(step_len), Photon_Path(4,end)];
        end

        t_scatter_stat = t_scatter_stat(1:cnt,:);
        batch_scatter_stat = [batch_scatter_stat; t_scatter_stat];
    end

    % keep EXACTLY nExamplePaths per batch
    batch_scatter_stat = batch_scatter_stat(1:nExamplePaths,:);

    % append to global storage
    scatter_stat = [scatter_stat; batch_scatter_stat];

    disp(['Batch ', num2str(batch), ...
          ' collected ', num2str(size(batch_scatter_stat,1)), ' photons'])
end
 
%% Assign outputs
x_in = scatter_stat(:,1); y_in = scatter_stat(:,2); z_in = scatter_stat(:,3);
x_ot = scatter_stat(:,4); y_ot = scatter_stat(:,5); z_ot = scatter_stat(:,6);
s    = scatter_stat(:,7); w    = scatter_stat(:,8);
no_of_photons = size(scatter_stat,1);
M_raw = t_model.G.M_raw;

%% Save results (parameter-safe)
filename = sprintf(['ua1_%d_us1_%d_n1_%0.2f_' ...
                    'ua2_%d_us2_%d_n2_%0.2f_' ...
                    'ua3_%d_us3_%d_n3_%0.2f.mat'], ...
                    ua1, us1, n1, ua2, us2, n2, ua3, us3, n3);

sim_info.ua1 = ua1; sim_info.us1 = us1; sim_info.n1 = n1;
sim_info.ua2 = ua2; sim_info.us2 = us2; sim_info.n2 = n2;
sim_info.ua3 = ua3; sim_info.us3 = us3; sim_info.n3 = n3;
sim_info.nBatches = nBatches;
sim_info.nPhotonsPerBatch = nPhotonsReq;
sim_info.nExamplePaths = nExamplePaths;

save(filename, 'scatter_stat', 't_model', 'M_raw', 'sim_info');

end

%% ------------------------------
%% Geometry function
function M = geometryDefinition(X,Y,Z,~)
    x_cyl_min = -7.5; x_cyl_max = 7.5;
    y_cyl = 0; z_cyl = 1.0; r_cyl = 0.3;

    M = ones(size(X));      % Layer 1 (dye)
    M(Z > 0.5) = 2;         % Layer 2 (tissue)
    
    mask = (X >= x_cyl_min) & (X <= x_cyl_max) & ...
           ((Y - y_cyl).^2 + (Z - z_cyl).^2 <= r_cyl^2);
    M(mask) = 3;            % Layer 3 (ink)
end

%% ------------------------------
%% Media properties function
function mediaProperties = mediaPropertiesFunc(var)
    mediaProperties = MCmatlab.mediumProperties;

    j = 1;
    mediaProperties(j).name = 'Layer 1 (dye)';
    mediaProperties(j).mua  = var{1};
    mediaProperties(j).mus  = var{2};
    mediaProperties(j).g    = 0.9;
    mediaProperties(j).n    = var{3};

    j = 2;
    mediaProperties(j).name = 'Layer 2 (tissue)';
    mediaProperties(j).mua  = var{4};
    mediaProperties(j).mus  = var{5};
    mediaProperties(j).g    = 0.9;
    mediaProperties(j).n    = var{6};

    j = 3;
    mediaProperties(j).name = 'Layer 3 (ink)';
    mediaProperties(j).mua  = var{7};
    mediaProperties(j).mus  = var{8};
    mediaProperties(j).g    = 0.9;
    mediaProperties(j).n    = var{9};
end
