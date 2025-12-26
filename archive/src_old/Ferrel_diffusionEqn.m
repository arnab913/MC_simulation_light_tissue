%% Farrell dipole reflectance for 3 cases: no dye, dye1, dye2
% Assumes LayerPath{1}, LayerPath{6}, LayerPath{11} already exist

%clear;  % keep in mind this will remove LayerPath if run standalone

% If LayerPath is already in workspace, comment out the clear above.

% --------------------------------------------------
% 1) Effective musp from Li for k = 1, 6, 11
% --------------------------------------------------
ks = [1 6 11];  % no dye, dye1, dye2
mu_sp_eff_all = zeros(size(ks));

for idx = 1:numel(ks)
    k  = ks(idx);
    Li = LayerPath{k}.Li;   % assume length >= 4
    L1 = Li(1);
    L2 = Li(2);
    L3 = Li(3);
    L4 = Li(4);

    switch k
        case 1   % no dye: all layers musp = 30 cm^-1
            mu_sp_eff_all(idx) = (L1 + L2 + L3 + L4) / ...
                                 (L1/30 + L2/30 + L3/30 + L4/30);

        case 6   % dye1: example musp: [5,10,20,30] cm^-1
            mu_sp_eff_all(idx) = (L1 + L2 + L3 + L4) / ...
                                 (L1/5 + L2/10 + L3/20 + L4/30);

        case 11  % dye2: example musp: [2,5,10,30] cm^-1
            mu_sp_eff_all(idx) = (L1 + L2 + L3 + L4) / ...
                                 (L1/2 + L2/5 + L3/10 + L4/30);

        otherwise
            error('Unexpected k = %d', k);
    end
end

mu_sp_eff_no_dye = mu_sp_eff_all(1);  % for k = 1
mu_sp_eff_dye1   = mu_sp_eff_all(2);  % for k = 6
mu_sp_eff_dye2   = mu_sp_eff_all(3);  % for k = 11

% --------------------------------------------------
% 2) Farrell dipole parameters
% --------------------------------------------------
mua = 0.05;               % uniform absorption [1/cm]
n   = 1.4;                % refractive index
r   = linspace(0.01, 6, 400);   % radial distance [cm], avoid r=0

% --------------------------------------------------
% 3) Compute reflectance for each case
% --------------------------------------------------
R_no_dye = dipole_reflectance(r, mua, mu_sp_eff_no_dye, n);
R_dye1   = dipole_reflectance(r, mua, mu_sp_eff_dye1,   n);
R_dye2   = dipole_reflectance(r, mua, mu_sp_eff_dye2,   n);

% Normalize each curve to its own maximum (optional)
R_no_dye = R_no_dye / max(R_no_dye);
R_dye1   = R_dye1   / max(R_dye1);
R_dye2   = R_dye2   / max(R_dye2);

% --------------------------------------------------
% 4) Plot: Farrell distribution of reflectance for all 3 cases
% --------------------------------------------------
figure; hold on;
semilogy(r, R_no_dye, 'k-', 'LineWidth', 2, 'DisplayName','No dye');
semilogy(r, R_dye1,   'r-', 'LineWidth', 2, 'DisplayName','Dye 1');
semilogy(r, R_dye2,   'b-', 'LineWidth', 2, 'DisplayName','Dye 2');

xlabel('Radial distance r [cm]');
ylabel('Diffuse reflectance (normalized)');
title('Farrell diffusion dipole reflectance (L_i-based musp_{eff})');
grid on;
ylim([1e-6 1]);
xlim([0 6]);
legend('show');
set(gca,'FontSize',12);

% --------------------------------------------------
% 5) Dipole reflectance function (Farrell et al.)
% --------------------------------------------------
function R = dipole_reflectance(r, mua, musp, n)
    % Steady-state diffusion dipole reflectance
    % r in cm, mua and musp (reduced scattering) in cm^-1, n dimensionless

    D     = 1 ./ (3 * (mua + musp));      % diffusion coefficient
    mueff = sqrt(mua ./ D);               % effective attenuation

    % Effective reflection coefficient Reff (Farrell approximation)
    Reff = -1.440./n.^2 + 0.710./n + 0.668 + 0.0636.*n;
    A    = (1 + Reff) ./ (1 - Reff);

    z0 = 1 ./ (mua + musp);              % source depth
    zb = 2 .* A .* D;                    % extrapolated boundary

    r1 = sqrt(r.^2 + z0.^2);
    r2 = sqrt(r.^2 + (z0 + 2*zb).^2);

    R = (1/(4*pi)) .* ( ...
          z0 .* (mueff + 1./r1) .* exp(-mueff .* r1) ./ r1.^2 + ...
         (z0 + 2*zb) .* (mueff + 1./r2) .* exp(-mueff .* r2) ./ r2.^2 );
end
