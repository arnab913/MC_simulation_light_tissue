function [x_in, y_in, z_in, x_ot, y_ot, z_ot, s, w, no_of_photons, M_raw] = Dyetissue_simulation(opt_bckg,cmd_size, opt_sgnl, ...
                  zSurface,wvlngth,nPhotonsReq,nExamplePaths, ...
                  beam_X,beam_Y,beam_phi,beam_tht)
%% Geometry definition
model = MCmatlab.model;
dl = 0.05; % spatial resolution (cm)
model.G.nx                = round(cmd_size(1)./dl);  % Number of bins in the x direction
model.G.ny                = round(cmd_size(2)./dl);  % Number of bins in the y direction
model.G.nz                = round(cmd_size(3)./dl);  % Number of bins in the z direction
model.G.Lx                = cmd_size(1); % [cm] x size of simulation cuboid
model.G.Ly                = cmd_size(2); % [cm] y size of simulation cuboid
model.G.Lz                = cmd_size(3); % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;  % Media properties defined as a function at the end of this file
model.G.mediaPropParams     = {opt_bckg(1),opt_bckg(2),opt_bckg(3),opt_bckg(4),opt_sgnl(1),opt_sgnl(2),opt_sgnl(3),opt_sgnl(4)}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc ,opt_sgnl(1),opt_sgnl(2),opt_sgnl(3),opt_sgnl(4)
model.G.geomFunc            = @geometryDefinition;   % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.
%model.G.geomFuncParams      = {zSurface,2,3}; % A cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc ,sig_pos(1),sig_pos(2),sig_pos(3),sig_size(1),sig_size(2),sig_size(3)
% model = plot(model,'G');
clearvars opt_bckg opt_sgnl cmd_size dl   %opt_sgnl sig_pos sig_size

%% Monte Carlo simulation
model.MC.simulationTimeRequested  = 1;       % [min] Time duration of the simulation
model.MC.matchedInterfaces        = false;   % Assumes all refractive indices are the same
model.MC.boundaryType             = 1;       % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.nPhotonsRequested = nPhotonsReq;    % The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.
model.MC.nExamplePaths = 5000;        % The code stores the paths of the first N photons for subsequent visualization during the plotting.
model.MC.wavelength = wvlngth;               % [nm] Excitation wavelength, used for determination of optical properties for excitation light
clearvars wvlngth nPhotonsReq

%% For a pencil beam, the "focus" is just a point that the beam goes through, here set to be the center of the cuboid:
model.MC.lightSource.sourceType   = 0;    % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.silentMode = false;
model.MC.lightSource.xFocus       = beam_X;   % [cm] x position of focus
model.MC.lightSource.yFocus       = beam_Y;   % [cm] y position of focus
model.MC.lightSource.zFocus       = zSurface; % [cm] z position of focus
model.MC.lightSource.theta        = beam_tht; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = beam_phi; % [rad] Azimuthal angle of beam center axis
clearvars zSurface beam_phi beam_tht beam_X beam_Y
scatter_stat = [];
no_of_photons = 0;

while size(scatter_stat,1) < nExamplePaths

    t_model = runMonteCarlo(model);
    % DEBUG: list fields in t_model.MC
    %disp(fieldnames(t_model.MC));

    %if isfield(t_model.MC, 'examplePaths')
    paths = t_model.MC.examplePaths;
    % else
    % error('MC.examplePaths not present. Check model.MC.nExamplePaths and MCmatlab version.');
    % end
   % paths = t_model.MC.examplePaths;

    % Find NaN separators
    spratrs = find(isnan(paths(1,:)));

    % Preallocate (max possible photons in this batch)
    t_scatter_stat = nan(length(spratrs)-1, 8);
    cnt = 0;   % counter for valid photons

    % ---- CLEAN LOGIC LOOP ----
    for i = 1:length(spratrs)-1

        start_col = spratrs(i)   + 1;
        end_col   = spratrs(i+1) - 1;

        % Skip empty photons
        if start_col > end_col
            continue
        end

        Photon_Path = paths(:, start_col:end_col);

        cnt = cnt + 1;

        % Compute statistics
        step_vecs = Photon_Path(1:3,2:end) - Photon_Path(1:3,1:end-1);
        step_len  = sqrt(sum(step_vecs.^2,1));

        t_scatter_stat(cnt,:) = [ ...
            Photon_Path(1,1), ...
            Photon_Path(2,1), ...
            Photon_Path(3,1), ...
            Photon_Path(1,end), ...
            Photon_Path(2,end), ...
            Photon_Path(3,end), ...
            sum(step_len), ...
            Photon_Path(4,end) ...
        ];
    end

    % Keep only filled rows
    t_scatter_stat = t_scatter_stat(1:cnt,:);

    scatter_stat = [scatter_stat; t_scatter_stat];

    disp('****************************************************')
    disp(['progress ... ', num2str(size(scatter_stat,1)/nExamplePaths*100), '%'])
    disp('****************************************************')

end

% %% These lines will run the Monte Carlo simulation with the provided parameters and subsequently plot the results:
% % figure(17)
% scatter_stat = [];
% no_of_photons = 0;
% while size(scatter_stat,1) < nExamplePaths
%     t_model = runMonteCarlo(model);
%     M_raw = t_model.G.M_raw;
%     % t_model = plot(t_model,'MC');
%     % do calc
%     spratrs = find(isnan(t_model.MC.examplePaths(1,:)));
%     srt_pnts =spratrs(i)+1;
%     fnsh_pnts =spratrs (i+1)-1;
%     strt_pnts = spratrs(2:2:end); strt_pnts = strt_pnts(1:end-1);
%     fnsh_pnts = spratrs(1:2:end); fnsh_pnts = fnsh_pnts(2:end-0);
%     t_scatter_stat = nan(length(fnsh_pnts),3+3+1+1); clearvars spratrs
%     for idx = 1:length(strt_pnts)
%         Photon_Path = t_model.MC.examplePaths(:,strt_pnts(idx)+1:fnsh_pnts(idx)-1);
%         % plot3(Photon_Path(1,:),Photon_Path(2,:),-Photon_Path(3,:)), hold on
%         t_scatter_stat(idx,:) = [...
%             (Photon_Path(1,1)), ...
%             (Photon_Path(2,1)), ...
%             (Photon_Path(3,1)), ...
%             Photon_Path(1,end), ...
%             Photon_Path(2,end), ...
%             Photon_Path(3,end), ...
%             sum(sqrt(sum((Photon_Path(1:3,2:end-0)-Photon_Path(1:3,1:end-1)).^2,1))), ...
%             Photon_Path(4,end), ...
%             ];        
%         clearvars Photon_Path
%     end
%     clearvars strt_pnts fnsh_pnts idx t_model
% 
%     no_of_photons = no_of_photons + size(t_scatter_stat,1);
%     scatter_stat = [scatter_stat ; t_scatter_stat];
%     disp('****************************************************')
%     disp(['progress ... ',num2str(size(scatter_stat,1)/nExamplePaths*100),'%'])
%     disp('****************************************************')
%     clearvars surf_photons t_scatter_stat
% end
% % figure(17), hold off, axis equal, view([-90 0]), grid on
% % axis([min(model.G.x) max(model.G.x) min(model.G.y) max(model.G.y) -max(model.G.z) -min(model.G.z)])
% % xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',24)
% % figure(17), hold off, axis equal, view([-0 +0]), grid on
% % axis([-7.5 +7.5 -7.5 +7.5 -max(model.G.z) -min(model.G.z)])
% % xlabel('x (cm)'), ylabel('y (cm)'), zlabel('z (cm)'), set(gca,'fontsize',24)
% 
% % histcount
% x_in = scatter_stat(:,1);                                                 % x -> input
% y_in = scatter_stat(:,2);                                                 % y -> input
% z_in = scatter_stat(:,3);                                                 % y -> input
% x_ot = scatter_stat(:,4);                                                 % x -> output
% y_ot = scatter_stat(:,5);                                                 % y -> output
% z_ot = scatter_stat(:,6);                                                 % y -> output
% s = scatter_stat(:,7);                                                    % s -> true distance
% w = scatter_stat(:,8);                                                    % w -> intensity

% ---- ASSIGN FUNCTION OUTPUTS ----
x_in = scatter_stat(:,1);
y_in = scatter_stat(:,2);
z_in = scatter_stat(:,3);

x_ot = scatter_stat(:,4);
y_ot = scatter_stat(:,5);
z_ot = scatter_stat(:,6);

s = scatter_stat(:,7);
w = scatter_stat(:,8);

no_of_photons = size(scatter_stat,1);
M_raw = t_model.G.M_raw;   % from last MC run
%model_out = t_model;          % return the full model object
save('Ink_100cm_tissueOnly.mat', 't_model', 'scatter_stat', 'M_raw');
end

function M = geometryDefinition (X,Y,Z,parameters)
% Geometry function(s) (see readme for details)
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
    mask = (X >= x_cyl_min) & (X <= x_cyl_max) & ...
           ((Y - y_cyl).^2 + (Z - z_cyl).^2 <= r_cyl^2);

    M(mask) = 2;   % Media 2 is the ink

end

function mediaProperties = mediaPropertiesFunc (var)
mediaProperties = MCmatlab.mediumProperties;
% Put in your own media property definitions below at 800 nm
j = 1;
mediaProperties(j).name  = 'tissue';
mediaProperties(j).mua   = double(var{1}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{2}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{3}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{4}); % The refractive index

j = 2;
mediaProperties(j).name  = 'ink signal';   % mimicking tumor
mediaProperties(j).mua   = double(var{5}); % Absorption coefficient [cm^-1]
mediaProperties(j).mus   = double(var{6}); % Scattering coefficient [cm^-1]
mediaProperties(j).g     = double(var{7}); % Henyey-Greenstein scattering anisotropy
mediaProperties(j).n     = double(var{8}); % The refractive index
end
