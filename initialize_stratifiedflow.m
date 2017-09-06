function [ gal, bilin ] = initialize_stratifiedflow(param, m, pol)
%INITIALIZE_STRATIFIEDFLOW Constructs Galerkin object for the PDE:
%
%   d/dt u = E*Lap u + psi_z + J(psi, u)
%   d/dt Lap psi = E*Lap^2 psi - u_z + s_y + J(psi, Lap psi)
%   d/dt s = E/Pr*Lap s + B*psi_z - S*psi_y + J(psi, s)
%

a = @(k) 2*pi*k/param.gamma;

% rhs{i}{j} = equation i on the rhs, applied to variable j

rhs{1}{1} = @(z, u, du, d2u)...
            param.E*(d2u - a(m)^2*u);

rhs{1}{2} = @(z, p, dp, d2p, d3p, d4p)...
            dp - 1i*a(m)*param.B*p;

rhs{1}{3} = @(z, s, ds, d2s)...
            0;

rhs{2}{1} = @(z, u, du, d2u)...  
            -du;
        
rhs{2}{2} = @(z, p, dp, d2p, d3p, d4p) ...
            param.E*(d4p - 2*a(m)^2*d2p + a(m)^4*p);

            
rhs{2}{3} = @(z, s, ds, d2s)...  
            1i*a(m)*s;            

rhs{3}{1} = @(z, u, du, d2u)...
            0;        
        
rhs{3}{2} = @(z, p, dp, d2p, d3p, d4p)...
            param.B*dp - 1i*a(m)*param.S*p;

rhs{3}{3} = @(z, s, ds, d2s)...
            param.E/param.Pr*(d2s - a(m)^2*s);
        
% lhs{i}{j} = equation i on the lhs, applied to variable j

lhs{1}{1} = @(z, u) u;
lhs{1}{2} = @(z, p, dp, d2p) 0;
lhs{1}{3} = @(z, s) 0;

lhs{2}{1} = @(z, u) 0;
lhs{2}{2} = @(z, p, dp, d2p) d2p - a(m)^2*p ;
lhs{2}{3} = @(z, s) 0;

lhs{3}{1} = @(z, u) 0;
lhs{3}{2} = @(z, p, dp, d2p) 0;
lhs{3}{3} = @(z, s) s;

% Boundary conditions
% bc_l{i}{j} jth boundary condition for variable i at the left endpoint
% j=1 if 2nd order eqn for u_i, j=1,2 if 4th order equation for u_i
bc_l{1}{1} = @(u, du) u;

bc_l{2}{1} = @(p, dp, d2p, d3p) p;
bc_l{2}{2} = @(p, dp, d2p, d3p) dp;

bc_l{3}{1} = @(s, ds) s;

% same at the right endpoint
bc_r{1}{1} = @(u, du) du;

bc_r{2}{1} = @(p, dp, d2p, d3p) p;
bc_r{2}{2} = @(p, dp, d2p, d3p) d2p;

bc_r{3}{1} = @(s, ds) ds + param.Bi*s;

gal = Galerkin(pol, bc_l, bc_r, rhs, lhs, 0, 1);

% The bilinear part
bilin = @(m1,m2,x, u,du,d2u,d3u, v,dv,d2v,d3v)...
    {1i*a(m1)*u{2}.*dv{1} - 1i*a(m2)*du{2}.*v{1};...
     1i*a(m1)*u{2}.*(d3v{2} - a(m2)^2*dv{2}) - 1i*a(m2)*du{2}.*(d2v{2} - a(m2)^2*dv{2});... 
     1i*a(m1)*u{2}.*dv{3} - 1i*a(m2)*du{2}.*v{3}};

end

