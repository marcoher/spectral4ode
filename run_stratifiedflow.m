%% 1) Setup polynomial basis
clear all

N = 20;

tol = 1e-12;

pol = Legendre(N, tol);

%% 2) Find eigenfunctions
close all

% Parameters
param = struct('E', 0.01,...
               'B', 1,...
               'S', 0.1,...
               'Pr', 1,...
               'gamma', 1,...
               'Bi', 0);
           
names = {'u', '\psi', '\sigma'};
           
m = 1;

% Equations
gal = initialize_stratifiedflow(param, m, pol);

% Eigenfunctions
nev = 1;

TransitionUtils.plot_eigfcns(gal, nev, param.gamma, m, 'flat', names);

%% 3) Look at f_m(lamda) = max_k Re beta_{m,k}(lambda)
figure(1)
clf

param.S = 0;
param.E = 0.001;

m_= 1:4;
E_ = 2*1e-3:1e-3:2*1e-2;
B_ = 0:0.1:4;
S_ = [0, 0.1, 1];%, 5, 10];

beta_ = zeros(numel(m_), numel(B_));
legendInfo = cell(numel(m_), 1);

for i=1:numel(m_)
    for j=1:numel(B_)
        param.B = B_(j);
        gal = initialize_stratifiedflow(param, m_(i), pol);
        beta_(i,j) = gal.eigenvalues(1, 'descend');
    end
    legendInfo{i} = ['m = ', num2str(m_(i))];
    plot(B_, real(beta_(i,:)));
    hold on
end

legend(legendInfo);
xlabel('$B$');
ylabel('$\max_k$ Re $\beta_{m,k}$', 'FontSize',12);
grid on

%% 4) Find critical lambda

clc
m = 2;
iters = 100;

for i=1:numel(S_)
    param.S = S_(i);
    B_c = TransitionUtils.critical_lambda(@initialize_stratifiedflow, pol, m,...
                                          param, 'B', 0, 10, iters);
    fprintf('Critical B = %0.6f \n', B_c);
end

%% 5) Neutral stability curves: compute
tic
iters = 100; % change to 500
lambdaTol = 1e-4; % change to 1e-10
betaTol = 1e-4; % change to 1e-10

B_c_m = zeros(numel(S_), numel(m_), numel(E_));
B_c = zeros(numel(S_), numel(E_));
m_c = zeros(numel(S_), numel(E_)); 

for i=1:numel(S_)
    param.S = S_(i);
    fprintf(['Working on S = ', num2str(S_(i)), '...\n']);
    [B_c_m(i,:,:), B_c(i,:), m_c(i,:), exit_flags] = TransitionUtils.neutral_stability_curves(@initialize_stratifiedflow, pol, m_,...
                                                        param, 'E', E_,...  
                                                        'B', 0, 50, 'false', ...
                                                        iters, lambdaTol, betaTol);
    fprintf(['Finished S = ', num2str(S_(i)), '.\n\n']);                                                
end

wtime = toc;
fprintf ('Computation completed in  %f seconds.\n', wtime );
%% 6) Neutral stability curves: plot
figure(2)
clf

for i=1:numel(S_)
    plot(E_, B_c(i,:), 'DisplayName', ['$S=', num2str(S_(i)), '$'])
    hold on
end
xlabel('$E$');
ylabel('$B^*$');
l1 = legend('show');
set(l1, 'interpreter', 'latex');


figure(3)
clf

i = 1;
legendInfo = cell(numel(m_),1);
for j=1:numel(m_)
    plot(E_, squeeze(B_c_m(i,j,:)))
    hold on
    legendInfo{j} = ['$m=',num2str(m_(j)),'$'];
end
xlabel('$E$');
ylabel('$B_m^*$');
%axis([min(alpha_) max(alpha_) 20 100]);
legend(legendInfo, 'interpreter', 'latex');

%% 7) Bifurcation numbers: compute

eta_ = zeros(numel(S_), numel(E_));
asymp = zeros(numel(S_), numel(E_));

for i=1:numel(S_)
    param.S = S_(i);
    for j=1:numel(E_)
        B = B_c(i,j);
        m = m_c(i,j);
        if j<size(m_c,2)
            if m_c(i,j+1)~=m_c(i,j)
                asymp(i,j) = 1;
            end
        end
        param.B = B;% = setfield(param, 'Re', Re);
        param.E = E_(j);% = setfield(par, 'b', b_(j));
        [gal_m, bilin] = initialize_stratifiedflow(param, m, pol);
        gal_0 = initialize_stratifiedflow(param, 0, pol);
        gal_2m = initialize_stratifiedflow(param, 2*m, pol);
        [eta, beta, u, h11, h20] = TransitionUtils.transition(gal_m, gal_0, gal_2m, bilin, m);
        eta_(i,j) = eta;
    end
    fprintf(['Found bifurcation numbers for S = ', num2str(S_(i)), '.\n']); 
end    

%% 8) Bifurcation numbers: plot
figure(4)
clf
figure(5)
clf

h = zeros(numel(S_),1);
h2 = zeros(numel(S_),1);

for i=1:numel(S_)
    re_eta = real(eta_(i,:));
    im_eta = abs(imag(eta_(i,:)));

    asymp_ = find(asymp(i,:)==1);
    
    figure(4)
    h(i) = plot(E_, re_eta, 'DisplayName', ['$S_',num2str(i),'=$',num2str(S_(i))]);
    hold on
    
    figure(5)
    h2(i) = plot(E_, im_eta, 'DisplayName', ['$S_',num2str(i),'=$',num2str(S_(i))]);
    hold on
end

figure(4)
xlabel('$E$');
ylabel('Re $\eta$');
l3=legend(h);
set(l3, 'interpreter', 'latex');

figure(5)
xlabel('$E$');
ylabel('$|$Im $\eta|$');
l4=legend(h2);
set(l4, 'interpreter', 'latex');

%% 9) Show bifurcated solution
close all

m = m_c(10);
param.B = B_c(10)+1;
param.E = E_(10);
    
[gal_m, bilin] = initialize_stratifiedflow(param, m, pol);
gal_0 = initialize_stratifiedflow(param, 0, pol);
gal_2m = initialize_stratifiedflow(param, 2*m, pol);

[eta, beta, u, h11, h20] = TransitionUtils.transition(gal_m, gal_0, gal_2m, bilin, m);

% Basic profile
u0 = @(y,z) {-param.B*z, 0, -param.B*y - param.S/param.B*z};

dt = 0.1;

sigma = imag(beta) - imag(eta)/real(eta) * real(beta);
T = min(2*2*pi/abs(sigma), 100);

M = TransitionUtils.animate_bifurcated_soln(beta, eta, param.gamma, m, dt, T, u, gal_m, h20, gal_2m, h11, gal_0, 'flat', u0, names);

[h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions
hf = figure; 
set(hf, 'position', [150 150 w h]);
axis off

movie(hf, M);

mplay(M, 10)