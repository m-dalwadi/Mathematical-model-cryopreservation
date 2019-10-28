function [t,xi_1,xi_2,xi_3,x_c, y_c, x_l, y_l, T_lc, T_ll, T_s, r_c, r_f, r_b] = nonlinear_transform_solution_switch_BC(cool,T_max,T_init,lambda)
%nonlinear_transform_solution Solve the dimensionless problem using a
%nonlinear transformation for xi_2 and by switching the interfacial
%condition forcing as suggested by the asymptotic analysis.

% cool is the dimensionless cooling rate of the system (beta in the paper).

%T_init is the initial time of the simulation - this affects the initial
%conditions, which are important if there is a boundary layer in the
%system for early time.

%T_max is the end time of the simulation.

%lambda controls how concentrated the grid points are in the liquid region
%near the freezing front. The higher lambda is, the more points near the
%freezing front. As lambda -> 0, we tend towards a uniform grid.

% Note that as cool gets higher, we want to start T_init earlier to capture
% all of the interesting behaviour. This also requires us to increase
% lambda, for the same reason. However, there is a difficulty associated
% with this - the boundary layer for x_l is much smaller than for T_ll.
% This means that we need to be quite careful as cool gets higher. As cool
% gets lower, everything tends towards a spatially uniform state, so we
% only need a coarse grid.

[p] = parameters(cool);

[p_ND] = ND_parameters(p);

n_1 = 80;
n_2 = 300;
n_3 = 80;

% Total number of points in the system
N = 3*n_1 + 3*n_2 + n_3 + 2;

xi_1 = linspace(0,1,n_1);
xi_2 = linspace(0,1,n_2);
xi_3 = linspace(0,1,n_3);

options = odeset('RelTol',1e-10,'AbsTol',1e-10);

[Init] = IC_after_transform(p_ND,T_init,xi_1,xi_2,xi_3,lambda);

tic
[t,Y] = ode15s(@(t,Y) discretization_ode(t,Y,xi_1,xi_2,xi_3,n_1,n_2,n_3,p_ND,lambda),[T_init T_max],Init,options);
toc

x_c = Y(:,1:n_1);
y_c = Y(:,n_1 + 1:2*n_1);
x_l = Y(:,2*n_1 + 1:2*n_1 + n_2);
y_l = Y(:,2*n_1 + n_2 + 1:2*n_1 + 2*n_2);

T_lc = Y(:,2*n_1 + 2*n_2 + 1: 3*n_1 + 2*n_2);
T_ll = Y(:,3*n_1 + 2*n_2 + 1: 3*n_1 + 3*n_2);
T_s = Y(:,3*n_1 + 3*n_2 + 1: 3*n_1 + 3*n_2 + n_3);

r_c = Y(:,3*n_1 + 3*n_2 + n_3 + 1);
r_f = Y(:,3*n_1 + 3*n_2 + n_3 + 2);

r_b = ones(size(t));







function dY = discretization_ode(t,Y,xi_1,xi_2,xi_3,n_1,n_2,n_3,p_ND,lambda)

%% The grid size for each domain
h_1 = xi_1(2) - xi_1(1);
h_2 = xi_2(2) - xi_2(1);
h_3 = xi_3(2) - xi_3(1);

%% Define concentration variables
x_c = Y(1:n_1);
y_c = Y(n_1 + 1:2*n_1);
x_l = Y(2*n_1 + 1:2*n_1 + n_2);
y_l = Y(2*n_1 + n_2 + 1:2*n_1 + 2*n_2);

%% Define temperature variables
% Note that T_l is split into cell and water parts. This means that the
% last element of T_lc is equivalent to the first element of T_ll.
T_lc = Y(2*n_1 + 2*n_2 + 1: 3*n_1 + 2*n_2);
T_ll = Y(3*n_1 + 2*n_2 + 1: 3*n_1 + 3*n_2);
T_s = Y(3*n_1 + 3*n_2 + 1: 3*n_1 + 3*n_2 + n_3);

r_c = Y(3*n_1 + 3*n_2 + n_3 + 1);
r_f = Y(3*n_1 + 3*n_2 + n_3 + 2);

%% Define system radius
r_b = 1;

%% Defining dr_f/dt

% Define flux into interface thru liquid. This corresponds to D x_r +
% \dot{r_f} x = 0 on r = r_f. Asymptotic analysis suggests that this
% controls r_f.
X_l_flux = 3*x_l(end) - 4*x_l(end-1) + x_l(end-2);
X_l_flux = X_l_flux /2 /h_2 / (r_f - r_c);
X_l_flux = X_l_flux *(exp(lambda) - 1)/lambda;
% Transformation to T = g/r
X_l_flux = X_l_flux./x_l(end);

dr_f = p_ND.D_x_water*(1/r_f - X_l_flux);

%% Defining dr_c/dt
% This corresponds to \dot{r_c} = -\kappa [(x_l-x_c) + sigma(y_l-y_c)].
dr_c = p_ND.sigma_x*(x_l(1) - x_c(end)) + p_ND.sigma_y*(y_l(1) - y_c(end));
dr_c = -p_ND.kappa * (1 + p_ND.nu *T_lc(end)/r_c) * dr_c;
% Transformation to x = g/r
dr_c = dr_c/r_c;

%% Defining dx_c/dt

% dx_c/dr = 0 at r = 0. This is equivalent to g = 0 at r = 0 when making
% the transformation x = g/r.

dx_c = zeros(size(x_c));

% Governing equation for x_c.
for k = 2:n_1-1
    dx_c(k) = p_ND.D_x_cell*(x_c(k+1) - 2*x_c(k) + x_c(k-1))/((r_c)^2*h_1^2) + dr_c*xi_1(k)*(x_c(k+1) - x_c(k-1))./((r_c)*2*h_1);
end
% BC at r = r_c. This corresponds to D x_c + \dot{r_c} x_c = w*(x_l-x_c).
RHS_f_c = p_ND.w * (1 + p_ND.nu *T_lc(end)/r_c) * (x_l(1) - x_c(end)) /r_c;
x_c_plus = x_c(end-1) + moving_derivative(p_ND.D_x_cell,0,r_c,h_1,r_c,dr_c,x_c(end),RHS_f_c);

dx_c(n_1) = p_ND.D_x_cell*(x_c_plus - 2*x_c(n_1) + x_c(n_1-1))/((r_c)^2*h_1^2) + dr_c*(x_c_plus - x_c(n_1-1))./((r_c)*2*h_1);

%% Defining dy_c/dt

% dy_c/dr = 0 at r = 0. This is equivalent to g = 0 at r = 0 when making
% the transformation y = g/r.

dy_c = zeros(size(y_c));

% Governing equation for y_c.
for k = 2:n_1-1
    dy_c(k) = p_ND.D_y_cell*(y_c(k+1) - 2*y_c(k) + y_c(k-1))/((r_c)^2*h_1^2) + dr_c*xi_1(k)*(y_c(k+1) - y_c(k-1))./((r_c)*2*h_1);
end

% BC at r = r_c. This corresponds to D y_r + \dot{r_c} y = w_y*(y_l-y_c),
% but note that we take w_y = 0 for non-permeable ions.
% RHS_h_c = 0;
RHS_h_c = p_ND.wy * (1 + p_ND.nu *T_lc(end)/r_c) * (y_l(1) - y_c(end)) /r_c;
y_c_plus = y_c(end-1) + moving_derivative(p_ND.D_y_cell,0,r_c,h_1,r_c,dr_c,y_c(end),RHS_h_c);

dy_c(n_1) = p_ND.D_y_cell*(y_c_plus - 2*y_c(n_1) + y_c(n_1-1))/((r_c)^2*h_1^2) + dr_c*(y_c_plus - y_c(n_1-1))./((r_c)*2*h_1);


%%
% Some helpful definitions for non-uniform grid.

exp_lam = [exp(-lambda*h_2),exp(lambda*xi_2),exp(lambda*(1+h_2))];
e_l = (exp_lam(1:end-1) + exp_lam(2:end))/2;

%% Defining dT_s/dt
% We need to define \dot{T_s} before x_l due to the BC for x_l at r = r_f.

dT_s = zeros(size(T_s));

% BC at r = r_f. This corresponds to the Stefan condition, which controls
% the temperature flux from the ice, from our asymptotic analysis.
T_l_flux = r_f*(3*T_ll(end) - 4*T_ll(end-1) + T_ll(end-2))/(r_f - r_c)/2/h_2;
T_l_flux = T_l_flux*(exp(lambda) - 1)/lambda;
RHS_T_l = p_ND.S*r_f^2*dr_f + T_s(1) + p_ND.k*(T_l_flux - T_ll(end));
RHS_T_l = RHS_T_l*(1 - r_f)/r_f;
T_s_minus = T_s(2) - 2*h_3*RHS_T_l;

dT_s(1) = p_ND.k_s*(T_s(2) - 2*T_s(1) + T_s_minus)/((r_b - r_f)^2*h_3^2) + (dr_f)*(T_s(2) - T_s_minus)./((r_b - r_f)*2*h_3);

dr_b = 0;

% Governing equation
for k = 2:n_3-1
    dT_s(k) = p_ND.k_s*(T_s(k+1) - 2*T_s(k) + T_s(k-1))/((r_b - r_f)^2*h_3^2) + ((1 - xi_3(k))*dr_f + dr_b*xi_3(k))*(T_s(k+1) - T_s(k-1))./((r_b - r_f)*2*h_3);
end

% BC at r = r_b. First term here arises due to differentiating T_ll/r_f
% with respect to time.

if t < ((1 - p_ND.alpha)/p_ND.cool_rate)
    dT_s(end) = dr_b*T_s(end)/r_b - p_ND.cool_rate*r_b;
else 
    dT_s(end) = dr_b*T_s(end)/r_b;
end

%% Defining dx_l/dt

dx_l = zeros(size(x_l));

% BC at r = r_c. This corresponds to D (x_l)_r + \dot{r_c} x = w(x_l - x_c)
LHS_f_l = p_ND.w * (1 + p_ND.nu *T_lc(end)/r_c) * (x_l(1) - x_c(end)) /r_c;
x_l_minus = x_l(2) - moving_derivative_NL(p_ND.D_x_water,r_c,r_f,h_2,r_c,dr_c,x_l(1),LHS_f_l,lambda);

dx_l(1) = p_ND.D_x_water*exp_lam(2)*(1 - exp(-lambda))^2*(e_l(2)*x_l(2) - (e_l(1) + e_l(2))*x_l(1) + e_l(1)*x_l_minus)/((r_f - r_c)^2*lambda^2*h_2^2) + exp_lam(2)*((1./exp_lam(2) - exp(-lambda))*dr_c + dr_f*(1 - 1./exp_lam(2)))*(x_l(2) - x_l_minus)./(lambda*(r_f - r_c)*2*h_2);

% Governing equation
for k = 2:n_2-1
    dx_l(k) = p_ND.D_x_water*exp_lam(k+1)*(1 - exp(-lambda))^2*(e_l(k+1)*x_l(k+1) - (e_l(k) + e_l(k+1))*x_l(k) + e_l(k)*x_l(k-1))/((r_f - r_c)^2*lambda^2*h_2^2) + exp_lam(k+1)*((1./exp_lam(k+1) - exp(-lambda))*dr_c + dr_f*(1 - 1./exp_lam(k+1)))*(x_l(k+1) - x_l(k-1))./(lambda*(r_f - r_c)*2*h_2);
end

% BC at r = r_f. This corresponds to T_s = -alpha*x_l.
dx_l(n_2) = -dT_s(1)/p_ND.alpha;


%% Defining dy_l/dt

dy_l = zeros(size(y_l));

% BC at r = r_c. This corresponds to D (y)_r + \dot{r_c} y = w_y(y_l -
% y_c), noting that w_y = 0 for impermeable ions.
% LHS_h_l = 0;
LHS_h_l = p_ND.wy * (1 + p_ND.nu *T_lc(end)/r_c) * (y_l(1) - y_c(end)) /r_c;
y_l_minus = y_l(2) - moving_derivative_NL(p_ND.D_y_water,r_c,r_f,h_2,r_c,dr_c,y_l(1),LHS_h_l,lambda);

dy_l(1) = p_ND.D_y_water*exp_lam(2)*(1 - exp(-lambda))^2*(e_l(2)*y_l(2) - (e_l(1) + e_l(2))*y_l(1) + e_l(1)*y_l_minus)/((r_f - r_c)^2*lambda^2*h_2^2) + exp_lam(2)*((1./exp_lam(2) - exp(-lambda))*dr_c + dr_f*(1 - 1./exp_lam(2)))*(y_l(2) - y_l_minus)./(lambda*(r_f - r_c)*2*h_2);

% Governing equation
for k = 2:n_2-1
    dy_l(k) = p_ND.D_y_water*exp_lam(k+1)*(1 - exp(-lambda))^2*(e_l(k+1)*y_l(k+1) - (e_l(k) + e_l(k+1))*y_l(k) + e_l(k)*y_l(k-1))/((r_f - r_c)^2*lambda^2*h_2^2) + exp_lam(k+1)*((1./exp_lam(k+1) - exp(-lambda))*dr_c + dr_f*(1 - 1./exp_lam(k+1)))*(y_l(k+1) - y_l(k-1))./(lambda*(r_f - r_c)*2*h_2);
end

% BC at r = r_f. This corresponds to D (y_l)_r + \dot{r_c} y_l = 0.
RHS_h_l = 0;
y_l_plus = y_l(end-1) + exp(-lambda)*moving_derivative_NL(p_ND.D_y_water,r_c,r_f,h_2,r_f,dr_f,y_l(end),RHS_h_l,lambda);

dy_l(n_2) = p_ND.D_y_water*exp_lam(n_2+1)*(1 - exp(-lambda))^2*(e_l(n_2+1)*y_l_plus - (e_l(n_2) + e_l(n_2+1))*y_l(n_2) + e_l(n_2)*y_l(n_2-1))/((r_f - r_c)^2*lambda^2*h_2^2) + exp_lam(n_2+1)*((1./exp_lam(n_2+1) - exp(-lambda))*dr_c + dr_f*(1 - 1./exp_lam(n_2+1)))*(y_l_plus - y_l(n_2-1))./(lambda*(r_f - r_c)*2*h_2);



%% Defining dT_lc/dt

dT_lc = zeros(size(T_lc));

% T_lc = 0 at r = 0 (as T = g/r)

% Governing equation
for k = 2:n_1-1
    dT_lc(k) = p_ND.k_l*(T_lc(k+1) - 2*T_lc(k) + T_lc(k-1))/((r_c)^2*h_1^2) + (dr_c*xi_1(k))*(T_lc(k+1) - T_lc(k-1))./((r_c)*2*h_1);
end

% BC at r = r_c.
% dT_lc/dr = dT_ll/dr
T_lc_plus = T_lc(n_1-1) + ((1 - exp(-lambda))/lambda)*(r_c*h_1/h_2/(r_f-r_c))*(-3*T_ll(1) + 4*T_ll(2) - T_ll(3));
dT_lc(end) = p_ND.k_l*(T_lc_plus - 2*T_lc(end) + T_lc(end-1))/((r_c)^2*h_1^2) + (dr_c)*(T_lc_plus - T_lc(end-1))./((r_c)*2*h_1);


%% Defining dT_ll/dt

dT_ll = zeros(size(T_ll));

% BC at r = r_c
% T_lc = T_ll
dT_ll(1) = dT_lc(end);

% Governing equation
for k = 2:n_2-1
    dT_ll(k) = p_ND.k_l*exp_lam(k+1)*(1 - exp(-lambda))^2*(e_l(k+1)*T_ll(k+1) - (e_l(k) + e_l(k+1))*T_ll(k) + e_l(k)*T_ll(k-1))/((r_f - r_c)^2*lambda^2*h_2^2) + exp_lam(k+1)*((1./exp_lam(k+1) - exp(-lambda))*dr_c + dr_f*(1 - 1./exp_lam(k+1)))*(T_ll(k+1) - T_ll(k-1))./(lambda*(r_f - r_c)*2*h_2);
end

% BC at r = r_f. T_l = -alpha*x_l
dT_ll(end) = -p_ND.alpha*(dx_l(end));



dY = [dx_c;dy_c;dx_l;dy_l;dT_lc;dT_ll;dT_s;dr_c;dr_f];

end


    function [F] = moving_derivative(D,a,b,h,s,sdot,f,RHS)
        % When we have a flux at a moving boundary, we have a BC of the
        % form D dx/dr + x sdot = RHS. When we transform this from
        % spherical to Cartesian, and then from moving to fixed boundaries,
        % we obtain something like x_{N+1} = x_{N-1} +/- 2*h
        % *(b-a)*(STUFF)/D and this function calculates the second term on
        % the right.
        
        % D = diffusivity
        % a = LHS bdy
        % b = RHS bdy
        % h = spatial step
        % s = bdy on which we are evaluating BC
        % sdot = velocity of that bdy
        % f = dependent variable of interest evaluated at the boundary
        
            F = 2*h*(b-a)*(f*(D/s - sdot) + s*RHS)/D;

    end


    function [F] = moving_derivative_NL(D,a,b,h,s,sdot,f,RHS,lambda)
        % When we have a flux at a moving boundary, we have a BC of the
        % form D dx/dr + x sdot = RHS. When we transform this from
        % spherical to Cartesian, and then from moving to fixed boundaries,
        % we obtain something like x_{N+1} = x_{N-1} +/- 2*h
        % *(b-a)*(STUFF)/D and this function calculates the second term on
        % the right.
        
        % D = diffusivity
        % a = LHS bdy
        % b = RHS bdy
        % h = spatial step
        % s = bdy on which we are evaluating BC
        % sdot = velocity of that bdy
        % f = dependent variable of interest evaluated at the boundary
        
            F = 2*h*(b-a)*lambda*(f*(D/s - sdot) + s*RHS)/(D*(1 - exp(-lambda)));
            
    end






end

