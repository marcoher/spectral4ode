function [gal, bilin ] = initialize_heating(param, m, pol)
%initialize_equation Set up a system of ODEs
%   Detailed explanation goes here

a = @(k) 2*pi*k/param.gam;

% rhs{i}{j} = equation i on the rhs, applied to variable j
rhs{1}{1} = @(z, u, du, d2u, d3u, d4u)...
    d4u - 2*a(m)^2*d2u + a(m)^4*u ;

rhs{1}{2} = @(z, v, dv, d2v)...
    -1i*a(m).*param.R*param.H_fun(z).*v ;

rhs{2}{1} = @(z, u, du, d2u, d3u, d4u)...
    -1i*a(m)*param.R*param.N_fun(z).*u ;

rhs{2}{2} = @(z, v, dv, d2v)...
    d2v - a(m)^2*v ;

% lhs{i}{j} = equation i on the lhs, applied to variable j
lhs{1}{1} = @(z, u, du, d2u) d2u - a(m)^2*u ;
lhs{1}{2} = @(z, v) 0;
lhs{2}{1} = @(z, u, du, d2u) 0;
lhs{2}{2} = @(z, v) param.Pr*v;

% Boundary conditions
% bc_l{i}{j} jth boundary condition for variable i at the left endpoint
% j=1 if 2nd order eqn for u_i, j=1,2 if 4th order equation for u_i
bc_l{1}{1} = @(u, du, d2u, d3u) u;
bc_l{1}{2} = @(u, du, d2u, d3u) du;
bc_l{2}{1} = @(v, dv) dv;
% same at the right endpoint
bc_r{1}{1} = @(u, du, d2u, d3u) u;
bc_r{1}{2} = @(u, du, d2u, d3u) du;
bc_r{2}{1} = @(v, dv) v;

x_fun = @(x) (x+1)/2;

gal = Galerkin(pol, bc_l, bc_r, rhs, lhs, x_fun);

% The bilinear part
bilin=@(m1,m2,x, u,du,d2u,d3u, v,dv,d2v,d3v)...
    {1i*a(m1)*(d2u{1} - a(m1)^2*u{1}).*dv{1} - 1i*a(m2)*(d3u{1} - a(m1)^2*du{1}).*v{1};...
    param.Pr*(1i*a(m1)*u{1}.*dv{2} - 1i*a(m2)*du{1}.*v{2}) };

end