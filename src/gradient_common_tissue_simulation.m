function [x_in, y_in, z_in, x_ot, y_ot, z_ot, s, w, no_of_photons, M_raw, B_wds, B_ds, Aink_phys, Aink_simple, debugPathsCell] = ...
    gradient_common_tissue_simulation(us1, us2, us3, ua5, us5, n5)

%addpath(genpath('/home/kniknam/MCmatlab 4.4.9.0'))

%% =========================
%% Simulation parameters
nBatches       = 200;         % number of batches
nPhotonsReq    = 100000;      % photons per batch
nExamplePaths  = 5000;        % example paths per batch
cmd_size       = [15, 15, 5]; % simulation volume [cm]
zSurface       = 0;           % surface depth [cm]
wvlngth        = 760;         % wavelength [nm]
beam_X         = 0; beam_Y = 0; beam_phi = 0; beam_tht = 0;
dl             = 0.05;        % spatial resolution [cm]

%% =========================
%% Plan B outputs (voxelized; no need to save all paths)
computeBanana      = true;    % accumulate banana volumes during run
computeInkAbs      = true;    % accumulate absorbed energy in ink during run
saveDebugPaths     = true;    % save small subset of paths for sanity checks
nDebugPaths        = 100;     % number of full paths to save (not 5000)
storeDebugAsSingle = true;    % reduce debug storage

%% =========================
%% Initialize MC model
model = MCmatlab.model;
model = plot(model,'G');
model.G.nx = round(cmd_size(1)/dl);
model.G.ny = round(cmd_size(2)/dl);
model.G.nz = round(cmd_size(3)/dl);
model.G.Lx = cmd_size(1); model.G.Ly = cmd_size(2); model.G.Lz = cmd_size(3);

% Media properties
model.G.mediaPropertiesFunc = @mediaPropertiesFunc;
model.G.mediaPropParams = {us1, us2, us3, ua5, us5, n5};
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

%% =========================
%% Accumulators (Plan B)
if computeBanana
    B_wds = zeros(model.G.nx, model.G.ny, model.G.nz, 'single'); % sum(wseg*ds)
    B_ds  = zeros(model.G.nx, model.G.ny, model.G.nz, 'single'); % sum(ds)
else
    B_wds = [];
    B_ds  = [];
end

Aink_phys   = 0;  % preferred physical estimator
Aink_simple = 0;  % weight-drop estimator (sanity check)

if saveDebugPaths
    debugPathsCell = cell(nDebugPaths,1);
    debugCount = 0;
else
    debugPathsCell = {};
    debugCount = 0;
end

% Ink cylinder geometry MUST match geometryDefinition()
x_cyl_min = -7.5; x_cyl_max = 7.5;
y_cyl = 0; z_cyl = 1.0; r_cyl = 0.3;

% Units note: cmd_size, dl, X/Y/Z are in cm, so muaInk must be in cm^-1
muaInk = ua5;

%% =========================
%% Main loops (keep your structure)
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

            Photon_Path = paths(:, start_col:end_col);  % 4 x N
            if size(Photon_Path,2) < 2
                continue
            end

            cnt = cnt + 1;

            % segment vectors + lengths
            step_vecs = Photon_Path(1:3,2:end) - Photon_Path(1:3,1:end-1);
            step_len  = sqrt(sum(step_vecs.^2,1));
            Ltot      = sum(step_len);

            % store the summary stats exactly like you had
            t_scatter_stat(cnt,:) = [ ...
                Photon_Path(1,1), Photon_Path(2,1), Photon_Path(3,1), ...
                Photon_Path(1,end), Photon_Path(2,end), Photon_Path(3,end), ...
                Ltot, Photon_Path(4,end)];

            % ---------- Save a small subset of full paths (debug only) ----------
            if saveDebugPaths && debugCount < nDebugPaths
                debugCount = debugCount + 1;
                if storeDebugAsSingle
                    debugPathsCell{debugCount} = single(Photon_Path);
                else
                    debugPathsCell{debugCount} = Photon_Path;
                end
            end

            % ---------- Plan B: voxelized banana + absorbed-in-ink ----------
            r = Photon_Path(1:3,:);    % 3 x N (cm)
            wstep = Photon_Path(4,:);  % 1 x N

            for k = 1:(size(r,2)-1)
                ds = step_len(k);
                if ds <= 0, continue; end

                rm = 0.5*(r(:,k) + r(:,k+1));  % midpoint (cm)

                % Banana voxel index (assumes X,Y centered around 0; Z from 0..Lz)
                if computeBanana
                    ix = floor((rm(1) + model.G.Lx/2)/dl) + 1;
                    iy = floor((rm(2) + model.G.Ly/2)/dl) + 1;
                    iz = floor((rm(3))/dl) + 1;

                    if ix>=1 && ix<=model.G.nx && iy>=1 && iy<=model.G.ny && iz>=1 && iz<=model.G.nz
                        wseg = 0.5*(wstep(k) + wstep(k+1));
                        B_wds(ix,iy,iz) = B_wds(ix,iy,iz) + single(wseg * ds);
                        B_ds(ix,iy,iz)  = B_ds(ix,iy,iz)  + single(ds);
                    end
                end

                if computeInkAbs
                    inInk = (rm(1) >= x_cyl_min) && (rm(1) <= x_cyl_max) && ...
                            ((rm(2)-y_cyl)^2 + (rm(3)-z_cyl)^2 <= r_cyl^2);

                    if inInk
                        % simple weight-drop
                        dw = wstep(k) - wstep(k+1);
                        if dw > 0
                            Aink_simple = Aink_simple + dw;
                        end
                        % preferred physical estimator
                        Aink_phys = Aink_phys + wstep(k) * (1 - exp(-muaInk * ds));
                    end
                end
            end
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

%% =========================
%% Assign outputs
x_in = scatter_stat(:,1); y_in = scatter_stat(:,2); z_in = scatter_stat(:,3);
x_ot = scatter_stat(:,4); y_ot = scatter_stat(:,5); z_ot = scatter_stat(:,6);
s    = scatter_stat(:,7); w    = scatter_stat(:,8);
no_of_photons = size(scatter_stat,1);
M_raw = t_model.G.M_raw;

% normalize banana volumes for comparisons (optional but recommended)
if computeBanana
    if sum(B_wds(:)) > 0, B_wds = B_wds ./ sum(B_wds(:)); end
    if sum(B_ds(:))  > 0, B_ds  = B_ds  ./ sum(B_ds(:));  end
end

% trim debug paths
if saveDebugPaths
    debugPathsCell = debugPathsCell(1:debugCount);
end

%% =========================
%% Save results (FILENAME SAFE: keep decimals unique)
filename = sprintf('us1_%g_us2_%g_us3_%g_ua5_%g_us5_%g_n5_%g.mat', us1, us2, us3, ua5, us5, n5);
filename = regexprep(filename, '\.', 'p');  % make filesystem friendly: 0.05 -> 0p05

sim_info.ua1 = 0.05; sim_info.us1 = us1; sim_info.n1 = 1.4;
sim_info.ua2 = 0.05; sim_info.us2 = us2; sim_info.n2 = 1.4;
sim_info.ua3 = 0.05; sim_info.us3 = us3; sim_info.n3 = 1.4;
sim_info.ua4 = 0.05; sim_info.us4 = 30;  sim_info.n4 = 1.4;
sim_info.ua5 = ua5;  sim_info.us5 = us5; sim_info.n5 = n5;  % <-- fixed typo

sim_info.nBatches = nBatches;
sim_info.nPhotonsPerBatch = nPhotonsReq;
sim_info.nExamplePaths = nExamplePaths;

sim_info.grid.dl = dl;
sim_info.grid.cmd_size = cmd_size;
sim_info.grid.nx = model.G.nx; sim_info.grid.ny = model.G.ny; sim_info.grid.nz = model.G.nz;
sim_info.grid.Lx = model.G.Lx; sim_info.grid.Ly = model.G.Ly; sim_info.grid.Lz = model.G.Lz;

% store ink geometry to reproduce mask later
sim_info.ink.x_cyl_min = x_cyl_min;
sim_info.ink.x_cyl_max = x_cyl_max;
sim_info.ink.y_cyl = y_cyl;
sim_info.ink.z_cyl = z_cyl;
sim_info.ink.r_cyl = r_cyl;

save(filename, 'scatter_stat', 't_model', 'M_raw', 'sim_info', ...
              'B_wds','B_ds','Aink_phys','Aink_simple','debugPathsCell', '-v7.3');

end

%% ------------------------------
%% Geometry function
function M = geometryDefinition(X,Y,Z,~)
    x_cyl_min = -7.5; x_cyl_max = 7.5;
    y_cyl = 0; z_cyl = 1.0; r_cyl = 0.3;

    % Tissue everywhere by default
    M = ones(size(X));       % Layer 1
    M(Z > 0.1) = 2;          % Layer 2
    M(Z > 0.3) = 3;          % Layer 3
    M(Z > 0.5) = 4;          % Layer 4 (bulk tissue)

    % Ink cylinder (media 5)
    mask = (X >= x_cyl_min) & (X <= x_cyl_max) & ...
           ((Y - y_cyl).^2 + (Z - z_cyl).^2 <= r_cyl^2);

    M(mask) = 5;
end

%% ------------------------------
%% Media properties function
function mediaProperties = mediaPropertiesFunc(var)
    mediaProperties = MCmatlab.mediumProperties;

    j = 1;
    mediaProperties(j).name = 'Layer 1 (grad dye L1)';
    mediaProperties(j).mua  = 0.05;
    mediaProperties(j).mus  = var{1};
    mediaProperties(j).g    = 0.9;
    mediaProperties(j).n    = 1.4;

    j = 2;
    mediaProperties(j).name = 'Layer 2 (grad dye L2)';
    mediaProperties(j).mua  = 0.05;
    mediaProperties(j).mus  = var{2};
    mediaProperties(j).g    = 0.9;
    mediaProperties(j).n    = 1.4;

    j = 3;
    mediaProperties(j).name = 'Layer 3 (grad dye L3)';
    mediaProperties(j).mua  = 0.05;
    mediaProperties(j).mus  = var{3};
    mediaProperties(j).g    = 0.9;
    mediaProperties(j).n    = 1.4;

    j = 4;
    mediaProperties(j).name = 'Layer 4 (Tissue L4)';
    mediaProperties(j).mua  = 0.05;
    mediaProperties(j).mus  = 30;
    mediaProperties(j).g    = 0.9;
    mediaProperties(j).n    = 1.4;

    j = 5;
    mediaProperties(j).name = 'Ink (Media 5)';
    mediaProperties(j).mua  = var{4};
    mediaProperties(j).mus  = var{5};
    mediaProperties(j).g    = 0.9;
    mediaProperties(j).n    = var{6};
end
