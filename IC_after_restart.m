function [Init] = IC_after_restart(xi_2,xi_2_new, x, y, T, r)
%GENERATE_IC Generate initial conditions for the solution after the system
%is stopped the first time.


Init_r_c = r.r_c.i(end);
Init_r_f = r.r_f.i(end);

Init_x_c = x.x_c.i(end,:);
Init_y_c = y.y_c.i(end,:);
Init_T_lc = T.T_lc.i(end,:);
Init_T_s = T.T_s.i(end,:);

Init_x_l = spline(xi_2_new,x.x_l.i(end,:),xi_2);
Init_y_l = spline(xi_2_new,y.y_l.i(end,:),xi_2);
Init_T_ll = spline(xi_2_new,T.T_ll.i(end,:),xi_2);

Init = [Init_x_c, Init_y_c, Init_x_l, Init_y_l, Init_T_lc, Init_T_ll, Init_T_s, Init_r_c, Init_r_f];

end
