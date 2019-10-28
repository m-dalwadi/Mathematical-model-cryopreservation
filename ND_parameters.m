function [p_ND] = ND_parameters(p)
%ND_PARAMETERS The dimensionless variables in the full system. We use the
%time scaling t ~ p_s c_s R_b^2/k_s.

% [p] = parameters;

p_ND = [];

% Diffusivity of x in cell
p_ND.D_x_cell = p.D_x_cell * p.p_s * p.c_s / p.k_s;
% Diffusivity of y in cell
p_ND.D_y_cell = p.D_y_cell * p.p_s * p.c_s / p.k_s;

% Relection coefficient of X
p_ND.sigma_x = p.sigma_x;
% Reflection coefficient of Y
p_ND.sigma_y = p.sigma_y;

% Diffusivity of x in water
p_ND.D_x_water = p.D_x_water * p.p_s * p.c_s / p.k_s;
% Diffusivity of y in water
p_ND.D_y_water = p.D_y_water * p.p_s * p.c_s / p.k_s;

% Thermal diffusivity ratio
p_ND.k_l = p.k_l*p.p_s * p.c_s/(p.k_s*p.p_l * p.c_l);

% Thermal diffusivity ratio
p_ND.k_s = 1;

Temp_diff = p.T_f0 - p.T_final;

% Variation of permeability/conductance withe temperature.
p_ND.nu = 0;


% Dimensionless cooling rate
p_ND.cool_rate = p.cool * ( p.p_s * p.c_s / p.k_s) * p.R_b^2 / Temp_diff;

% Relative temperature difference due to solute
p_ND.alpha = p.alpha * p.X_0/(Temp_diff);

% Hydraulic conductivity to water to (?)
p_ND.kappa = p.kappa * ( p.p_s * p.c_s / p.k_s) * p.R_b * p.X_0 * p.R * p.T_f0;

% Membrane permeability to X to (?)
p_ND.w = p.w *  p.R * ( p.p_s * p.c_s / p.k_s) * p.R_b * p.T_f0;

% Membrane permeability to Y to (?)
p_ND.wy = p.wy *  p.R * ( p.p_s * p.c_s / p.k_s) * p.R_b * p.T_f0;

% Thermal diffusivity ratio
p_ND.k = p.k_l/p.k_s;

% Density ratio
p_ND.p = p.p_l/p.p_s;

% Stefan number
p_ND.S = p.L/Temp_diff/p.c_s;

% Initial cell radius
p_ND.gamma = p.R_c / p.R_b;

% Initial Y_e concentration
p_ND.Y_e = p.Y_e/p.X_0;
% Initial Y_i concentration
p_ND.Y_i = p.Y_i/p.X_0;



end

