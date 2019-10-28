function [t,xi, x, y, T, r,lambda,T_init,T_init_1,p,p_ND] = Full_solver(cool)
%Full_solver Solve the full cryopreservation system. For high cooling rates
% is starts off by solving using a boundary layer solver, then switches to
% a uniform grid.

% "cool" is the dimensional cooling rate, in K s^{-1}.

[p] = parameters(cool);
[p_ND] = ND_parameters(p);

T_max = 10/p_ND.w;


% If cooling rate is very quick, we start with a solver to resolve the
% initial boundary layer, then move to a uniform grid since a uniform grid
% ends up being much quicker to solve over the long time.
if cool > 5e-1
    
    T_max_i = 10;
    
    % lambda is a measure of how bunched up the points are towards the BL.
    % This step is to determine how large lambda needs to be for the given
    % cooling rate. We choose T_init such that cool*t/alpha << 1, and
    % lambda such that we have at least 10 points in the boundary layer. We
    % fix the fact that we use 300 grid points in the liquid region. 
    P = 3*sqrt(p_ND.alpha/p_ND.cool_rate)/10/(1-p_ND.gamma);
    
    ww = @(x) x./(exp(x) - 1) - P;
    
    l_0 = 1; % initial point
    lambda = fzero(ww,l_0);
    
    T_init = (lambda*(1 - p_ND.gamma)/3/(exp(lambda) - 1)).^2;
    
    [t.i,xi.xi_1.i,xi.xi_2.i,xi.xi_3.i,x.x_c.i, y.y_c.i, x.x_l.i, y.y_l.i, T.T_lc.i, T.T_ll.i, T.T_s.i, r.r_c.i, r.r_f.i, r.r_b.i] ...
        = nonlinear_transform_solution_switch_BC(cool,T_max_i,T_init,lambda);
    
    T_init_1 = T_max_i;
    
    xi_2_new = (1 - exp(-lambda*xi.xi_2.i))./(1 - exp(-lambda));
    n_2 = length(xi_2_new);

    [Init] = IC_after_restart(xi.xi_2.i,xi_2_new, x, y, T, r);
   
    [t.f,xi.xi_1.f,xi.xi_2.f,xi.xi_3.f,x.x_c.f, y.y_c.f, x.x_l.f, y.y_l.f, T.T_lc.f, T.T_ll.f, T.T_s.f, r.r_c.f, r.r_f.f, r.r_b.f] ...
        = nonlinear_transform_sol_after_BL(cool,n_2,T_max,T_init_1,Init);
    
    
elseif cool > 5e-4
    % We choose lambda to be very small for lower cooling rates, since this
    % corresponds to taking a uniform grid. As the cooling rate decreases,
    % the BL subsides and the simulation has to be started at later times,
    % otherwise the fact that the problem is almost spatially independent
    % means that the solver requires high accuracy to avoid numerical error
    % creeping in when evaluating derivatives.
    lambda = 1e-5;
    T_init = 1e-1;
    T_init_1 = T_max;
    [t,xi.xi_1,xi.xi_2,xi.xi_3,x.x_c, y.y_c, x.x_l, y.y_l, T.T_lc, T.T_ll, T.T_s, r.r_c, r.r_f, r.r_b] ...
        = nonlinear_transform_solution_switch_BC(cool,T_max,T_init,lambda);
    disp('No time split')
elseif cool > 5e-5
    % We choose lambda to be very small for lower cooling rates, since this
    % corresponds to taking a uniform grid. As the cooling rate decreases,
    % the BL subsides and the simulation has to be started at later times,
    % otherwise the fact that the problem is almost spatially independent
    % means that the solver requires high accuracy to avoid numerical error
    % creeping in when evaluating derivatives.
    lambda = 1e-5;
    T_init = 3e0;
    T_init_1 = T_max;
    [t,xi.xi_1,xi.xi_2,xi.xi_3,x.x_c, y.y_c, x.x_l, y.y_l, T.T_lc, T.T_ll, T.T_s, r.r_c, r.r_f, r.r_b] ...
        = nonlinear_transform_solution_switch_BC(cool,T_max,T_init,lambda);
    disp('No time split')
else
    lambda = 1e-5;
    % We choose lambda to be very small for lower cooling rates, since this
    % corresponds to taking a uniform grid. As the cooling rate decreases,
    % the BL subsides and the simulation has to be started at later times,
    % otherwise the fact that the problem is almost spatially independent
    % means that the solver requires high accuracy to avoid numerical error
    % creeping in when evaluating derivatives.
    T_init = 50;
    T_init_1 = T_max;
    [t,xi.xi_1,xi.xi_2,xi.xi_3,x.x_c, y.y_c, x.x_l, y.y_l, T.T_lc, T.T_ll, T.T_s, r.r_c, r.r_f, r.r_b] ...
        = nonlinear_transform_solution_switch_BC(cool,T_max,T_init,lambda);
    disp('No time split')
end


end

