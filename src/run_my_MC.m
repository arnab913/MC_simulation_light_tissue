% Optical properties: [mua, mus, g, n]
opt_bckg = [0.05, 30, 0.9, 1.4];   % background (muscle) tissue
opt_sgnl = [100, 30, 0.9, 1.4];   % same as muscle as no signal condition first.ua of ink is 1,7,30,100 cm-1

% Simulation volume [cm]
cmd_size = [10, 10, 5];

% Tumor center and size [cm]
% sig_pos  = [0, 0, 2];   % center at z = 2 cm
% sig_size = [1, 1, 1];   % 1 cm semi-axes

% Surface depth [cm]
zSurface = 0;

% Wavelength [nm]
wvlngth = 760;

% Photon counts
nPhotonsReq   = 1e5;
nExamplePaths = 5000; % change to nPhotonsReq later

% Beam position and orientation
beam_X  = 0;
beam_Y  = 0;
beam_phi = 0;
beam_tht = 0;

% Run simulation
[x_in,y_in,z_in, x_ot,y_ot,z_ot, s, w, no_of_photons, M_raw] = ...
    tissueOnly_simulation(opt_bckg,cmd_size, opt_sgnl, ...
                  zSurface,wvlngth,nPhotonsReq,nExamplePaths, ...
                  beam_X,beam_Y,beam_phi,beam_tht);


% opt_sgnl,sig_pos,sig_size variables are ommited as no signal inside