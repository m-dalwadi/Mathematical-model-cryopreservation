function [Init] = IC_after_transform(p_ND,T_init,xi_1,xi_2,xi_3,lambda)
%GENERATE_IC Generate initial conditions for the dimensionless solution
%using the asymptotic analysis in Appendix C.


Init_r_c = p_ND.gamma;

r_1 = xi_1.*Init_r_c;

A = 4*sqrt(p_ND.D_x_water)*p_ND.cool_rate/sqrt(pi)/p_ND.alpha/3;
Init_r_f = 1 - A*T_init.^(3/2);

xi_2_NL = (1 - exp(-lambda*xi_2))/(1 - exp(-lambda));
r_2 = xi_2_NL.*(Init_r_f - Init_r_c) + Init_r_c;

Init_r_b = 1;
r_3 = xi_3.*(Init_r_b - Init_r_f) + Init_r_f;

F = @(x) (2*x.^2 + 1).*erfc(x) - 2*x.*exp(-x.^2)./sqrt(pi);

Init_x_c = (r_1);
Init_y_c = p_ND.Y_i* (r_1);

rho_2 = (1 - xi_2)*(1 - Init_r_c)*lambda/(exp(lambda) - 1)/sqrt(4*T_init);

X = T_init*F(rho_2/sqrt(p_ND.D_x_water));
X = X*p_ND.cool_rate/p_ND.alpha;

Y = T_init*F(rho_2/sqrt(p_ND.D_y_water));
Y = Y*p_ND.cool_rate*sqrt(p_ND.D_x_water/p_ND.D_y_water)/p_ND.alpha;

Init_x_l = (r_2 + X);
Init_y_l = p_ND.Y_e*(r_2 + Y);


Init_T_lc =  r_1.*(-p_ND.alpha);

T_L = T_init.*F(rho_2/sqrt(p_ND.k_l));
T_L = p_ND.cool_rate*T_L;

Init_T_ll = r_2.*(-p_ND.alpha) - T_L;

Init_T_s = -p_ND.alpha - p_ND.cool_rate*T_init;
Init_T_s = Init_T_s.*r_3;

Init = [Init_x_c, Init_y_c, Init_x_l, Init_y_l, Init_T_lc, Init_T_ll, Init_T_s, Init_r_c, Init_r_f];



end
