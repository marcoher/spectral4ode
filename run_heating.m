%% 0) Setup polynomial basis
clear all
N = 50;
tol = 1e-10;
pol = Legendre(N, tol);

%% 1) Find a few eigenfunctions
close all

% Parameters
gam = 2;
delta = 1/100;
R = 53.4210;
Pr = 1;
H_fun = @(z) 1+delta*z;
%H_fun = @(z) 1;
N_fun=@(z) z;
%N_fun=@(z) (z+z.^2)/2;
%N_fun=@(z) (2*z+1/2*z.^3-3/2*z.^2);
%N_fun=@(z)z+(3-2*cos(2*pi*z)-cos(4*pi*z))/(4*pi);
%N_fun = @(z) (z-z.^2);
%N_fun = @(z) realpow(z, 0.2);
%N_fun = @(z) z + 0.5*sin(10*z.^2);
%N_fun = @(z) (1 - cos(20*pi*z))/2;
%N_fun = @(z) z.*exp(-5*(z-1));
%N_fun = @(z) exp(z) - z - 1;
%N_fun = @(z) 200*z.*exp(-10*z);
%N_fun = @(z) 1;

param = struct('R', R,...
               'Pr', Pr,...
               'gam', gam,...
               'H_fun', H_fun,...
               'N_fun', N_fun);
m = 1;
         
% Equations
gal = initialize_heating(param, m, pol);

% Eigenfunctions
nev = 1;
TransitionUtils.plot_eigfcns(gal, nev, gam, m, 'flat');


%% 2) Plot beta as function of lambda
close all

m_= 1:3;
R_ = 40:90;

beta_ = TransitionUtils.plot_beta_lambda(@initialize_heating, m_, param, 'R', R_, pol);

%% 3) Find critical lambda
close all
clc

m = 1;
iters = 5;

R_c = TransitionUtils.critical_lambda(@initialize_heating, pol, m,...
                                      param, 'R', 40, 90, iters);

fprintf('Critical R = %0.6f \n', R_c);

%% 4) Neutral stability curves
% these curves should look smooth - if they don't, increase/decrease
% iters/betaTol

close all

m_ = 1:4;
gam_ = 3:0.1:10;

iters = 4;
betaTol = 1e-4;

[R_c_, m_c_] = TransitionUtils.neutral_stability_curves(@initialize_heating, pol, m_,...
                                                        param, 'gam', gam_,...  
                                                        'R', 1, 500, ...
                                                        iters, betaTol);
xlabel('\gamma'); % looks better

%figure
%plot(gam_, R_c_) % if you just want the critical R for each gam

%% 5) Find transition number
close all
%clc
format long 

% These numbers are found beforehand
param.R = R_c;
m = 1;

[gal_m, bilin] = initialize_heating(param, m, pol);
gal_0 = initialize_heating(param, 0, pol);
gal_2m = initialize_heating(param, 2*m, pol);

[eta, beta, u, h11, h20] = TransitionUtils.transition(gal_m, gal_0, gal_2m, bilin, m);
fprintf('Transition number:\n %0.10f+i%0.10f \n', real(eta), imag(eta));

%% 6) Transition as function of a parameter
close all

iters = 10;
gam_ = 1:0.5:5;
bif_ = zeros(size(gam_));
for i=1:numel(gam_)
    param.gam = gam_(i);
    param.R = TransitionUtils.critical_lambda(@initialize_heating, pol, m, param, 'R', 10, 100, iters);
    
    [gal_m, bilin] = initialize_heating(param, m, pol);
    gal_0 = initialize_heating(param, 0, pol);
    gal_2m = initialize_heating(param, 2*m, pol);

    [eta, beta, u, h11, h20] = TransitionUtils.transition(gal_m, gal_0, gal_2m, bilin, m);
    bif_(i) = eta;
end

plot(gam_, bif_)

%% 7) Show bifurcated solution
close all

% Basic profile
u0 = @(z) {-1/(pi^4*param.E)*sin(pi*z)};

dt = 1;

sigma = imag(beta) - imag(eta)/real(eta) * real(beta);
T = min(2*pi/abs(sigma), 100);

M = TransitionUtils.animate_bifurcated_soln(beta, eta, 2/param.a, m, dt, T, u, gal_m, h20, gal_2m, h11, gal_0, 'flat', u0);
