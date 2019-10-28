function [p] = parameters(cool)
%Define the parameters for the cryopreservation problem.

p = [];

% Thermal parameters
p.k_l = 0.6; % Thermal conductivity of water (W K^(-1) m^(-1))
p.p_l = 1e3; % Density of water (kg m^(-3))
p.c_l = 4e3; % Specific heat capacity of water (J kg^(-1) K^(-1))

p.k_s = 2.2; % Thermal conductivity (W K^(-1) m^(-1))
p.p_s = 9e2; % Density (kg m^(-3))
p.c_s = 2e3; % Specific heat capacity (J kg^(-1) K^(-1))

p.L = 3.4e5; % Latent heat of freezing (J kg^(-1))
p.alpha = 4e-3; % (K m^3 mol^(-1))


% Diffusivities
p.D_x_water = 1e-9; % Diffusivity of x in water (m^2 s^(-1))
p.D_y_water = 2e-9; % Diffusivity of y in water (m^2 s^(-1))
p.D_x_cell = 2e-10; % Diffusivity of x in cell (m^2 s^(-1))
p.D_y_cell = 4e-10; % Diffusivity of y in cell (m^2 s^(-1))

% Cell membrane parameters
% p.gamma = 1e-4; % Surface tension of cell membrane (kg s^(-2))
p.kappa = 5e-15; % Hydraulic conductivity (m^2 s kg^(-1))
p.R_c = 5e-5; % Typical cell radius (m)
p.V_c = 4*pi*p.R_c.^3/3; % Typical cell volume (m^3)
p.sigma_x = 0.65; % Relection coefficient of X
p.sigma_y = 1;%0.6; % Reflection coefficient of Y
p.R = 8.3; % Universal gas coefficient (kg m^2 s^(-2) K^(-1) mol^(-1))
p.R_b = 5e-4; % System radius (m)
p.V_b = 4*pi*p.R_b.^3/3 - p.V_c; % System volume (not including cell volume) (m^3)
p.X_0 = 1e3; % Initial concentration of X (mol m^(-3))
% p.Y_e = 2.4e-11 / p.V_b; % Initial (and total) concentration of Y outside cell (mol m^(-3))
% p.Y_i = 4.1e-14 / p.V_c; % Initial (and total) concentration of Y inside cell (mol m^(-3))
p.Y_e = 1e2; % Initial (and total) concentration of Y outside cell (mol m^(-3))
p.Y_i = 1e2; % Initial (and total) concentration of Y inside cell (mol m^(-3))
% p.phi_x = 1; % Cell--water partition coefficient for X
p.w = 5e-14; % Permeability of membrane to X (s mol m^(-1) kg^(-1))
p.wy = 0;%5e-13; % Permeability of membrane to Y (s mol m^(-1) kg^(-1))
% p.T_0 = 277; % Initial temperature of system (4 degrees C)
p.T_f0 = 273; % Freezing temperature of water with no solute (0 degrees C)
p.T_final = 200; % Final temperature of the system
p.cool = cool; % Cooling rate of liquid (K s^(-1))




end

