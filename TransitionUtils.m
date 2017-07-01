classdef TransitionUtils
    %TransitionUtils Some methods for studying dynamic transitions using 
    %the Galerkin and Legendre/Chebyshev classes
    %   Detailed explanation goes here
    
    methods(Static)
        function [xx, zz, uu] = to_2D(u, a, x, gal, flag)
        %to_2D Undoes separation of variables
        %   returns uu(x,z) = exp( i*a*x )* u(z)    
        
            [xx, zz] = meshgrid(x, gal.x_fun(gal.pol.x));
            uu = cell(numel(u));
            for i=1:numel(u)
                uu{i} = zeros(size(xx));
                for j=1:numel(x)
                    uu{i}(:,j) = (cos(a*x(j))+1i*sin(a*x(j)))*u{i};
                end
            end
            if strcmp(flag, 'polar')
                for i=1:numel(uu)
                    [xx, zz, uu{i}] = pol2cart(xx, zz, uu{i});
                end
            end
        end
        
        function [bif, beta, u, h11, h20] = transition(gal_m, gal_0, gal_2m, bilin, m)
        %transition Finds bifurcation number and center manifold function
        %           for 2D separable equation

            % Solve the eigenvalue problem with m
            [beta, uu, duu, d2uu, d3uu] = gal_m.eigenfunctions(1, 'normal');

            if real(beta)<0
                warning('Mode %d is linearly stable.',m)
            end

            % Normalize in L^2 (i.e. sum_i int_{-1}^1 |u_i(x)|^2 dx = 1)
            norm = sqrt(gal_m.inner(uu{1}, uu{1}));
            
            % We only care about the first (normalized) eigenfunction
            u = cell(size(uu{1}));
            du = cell(size(duu{1}));
            d2u = cell(size(d2uu{1}));
            d3u = cell(size(d3uu{1}));
            
            v = cell(size(u));
            dv = cell(size(du));
            d2v = cell(size(d2u));
            d3v = cell(size(d3u));
            
            for i=1:numel(uu{1}) %faster than cellfun
                u{i} = uu{1}{i}/norm;
                v{i} = conj(u{i});
                du{i} = duu{1}{i}/norm;
                dv{i} = conj(du{i});
                d2u{i} = d2uu{1}{i}/norm;
                d2v{i} = conj(d2u{i});
                d3u{i} = d3uu{1}{i}/norm;
                d3v{i} = conj(d3u{i});
            end
            
            % Find G_(m,m) (U_m, U_m)
            xx = gal_m.x_fun(gal_m.pol.x);
            g20 = bilin(m, m, xx, u, du, d2u, d3u, u, du, d2u, d3u);

            % Find G_(m,-m) (U_m, U_{-m})+G_(-m,m)(U_{-m}, U_m)
            g11 = bilin(m, -m, xx, u, du, d2u, d3u, v, dv, d2v, d3v);
            g11_2 = bilin(-m, m, xx, v, dv, d2v, d3v, u, du, d2u, d3u);
            for n = 1:numel(g11)
                g11{n} = g11{n}+g11_2{n};
            end

            % Find phi_20
            [h20, dh20, d2h20, d3h20] = gal_2m.solve_linear(2*beta, g20);

            % Find phi_11
            [h11, dh11, d2h11, d3h11] = gal_0.solve_linear(2*real(beta), g11);

            % Compute G(U_m, phi_11)+G(phi_11,U_m)+G(conj U_m,phi_2m)+G(phi_2m,conj U_m)
            gfin = bilin(m, 0, xx, u, du, d2u, d3u, h11, dh11, d2h11, d3h11);
            gfin_2 = bilin(0, m, xx, h11, dh11, d2h11, d3h11, u, du, d2u, d3u);
            gfin_3 = bilin(2*m, -m, xx, h20, dh20, d2h20, d3h20, v, dv, d2v, d3v);
            gfin_4 = bilin(-m, 2*m, xx, v, dv, d2v, d3v, h20, dh20, d2h20, d3h20);
            for n=1:numel(gfin)
                gfin{n} = gfin{n}+gfin_2{n}+gfin_3{n}+gfin_4{n};
            end

            % Find dual eigenfunctions 
            [~,u_dual] = gal_m.eigenfunctions(1, 'dual');
            uu_dual = gal_m.BPhiT(u_dual);

            num = gal_m.inner(gfin, uu_dual);
            den = u_dual'*(gal_m.opLhs*gal_m.FPhiT(u))/gal_m.sc;
            bif = num/den;
        end
        
        function [] = plot_eigfcns(gal, nev, period, m, flag)
        %plot_eigfcns plots eigenfunctions
        
            [beta, u] = gal.eigenfunctions(nev, 'normal');
            x = (0:0.01:period)';
            for i=1:nev
                for j=1:gal.neq
                    figure(i)
                    subplot(2,gal.neq,j);
                    if strcmp(flag, 'flat')
                        gal.plot_as_flat(u{i}{j}, x, m, 'real');
                    elseif strcmp(flag, 'polar')
                        gal.plot_as_polar(u{i}{j}, x, m, 'real');
                    else
                        error('flag must be either ''flat'' or ''polar''.');
                    end
                    xlabel('x'); 
                    ylabel('z'); 
                    title(['\Re u^{(',num2str(j),')}_', num2str(i), ', \Re \beta_{', num2str(m), ',', num2str(i),'}=', num2str(beta(i))]);
                    
                    subplot(2,gal.neq,gal.neq+j);
                    if strcmp(flag, 'flat')
                        gal.plot_as_flat(u{i}{j}, x, m, 'imag');
                    elseif strcmp(flag, 'polar')
                        gal.plot_as_polar(u{i}{j}, x, m, 'imag');
                    else
                        error('flag must be either ''flat'' or ''polar''.');
                    end
                    xlabel('x'); 
                    ylabel('z'); 
                    title(['\Im u^{(',num2str(j),')}_', num2str(i), ', \Re \beta_{', num2str(m), ',', num2str(i),'}=', num2str(beta(i))]);
                    
                end
    
            end
        end
        
        function beta_ = plot_beta_lambda(initializer, m_, param, lambda_name, lambda_, pol)
        %plot_beta_lambda plots max_k Real beta_{m,k}(lambda)
        
            beta_ = zeros(numel(m_), numel(lambda_));
            legendInfo = cell(numel(m_), 1);
            for i=1:numel(m_)
                for j=1:numel(lambda_)
                    param = setfield(param, lambda_name, lambda_(j));
                    gal = initializer(param, m_(i), pol);
                    beta_(i,j) = gal.eigenvalues(1, 'descend');
                end
                legendInfo{i} = ['m = ', num2str(m_(i))];
                plot(lambda_, real(beta_(i,:)), 'color', rand(1,3));
                hold on
            end
            legend(legendInfo);
            xlabel(lambda_name);
            ylabel('max_k \Re \beta_{m,k}');
            grid on
        end
        
        function M = animate_bifurcated_soln(beta, eta, period, m, dt, T, u,...
                                             gal_m, h20, gal_2m, h11, gal_0, flag, u0_)
            rho = sqrt( -real(beta)/real(eta) );
            sigma = imag(beta) - imag(eta)/real(eta) * real(beta);
            am = 2*pi*m/period;
            xx = 0:0.01:period;
            
            [x_, z_, uu] = TransitionUtils.to_2D(u, am, xx, gal_m, flag);
            [~, ~, hh20] = TransitionUtils.to_2D(h20, 2*am, xx, gal_2m, flag);
            [~, ~, hh11] = TransitionUtils.to_2D(h11, 0, xx, gal_0, flag);

            u0 = u0_(z_);
            N = ceil(T/dt);
            M(N)= struct('cdata',[],'colormap',[]);
            fig = figure;
            for j=0:N-1
                z = rho*(cos(sigma*j*dt)+1i*sin(sigma*j*dt));
                for i=1:gal_m.neq
                    u_t = u0{i}+2*real(z*uu{i}+z^2*hh20{i})+abs(z)^2*hh11{i};
                    subplot(gal_m.neq, 1, i);
                    contourf(x_, z_, u_t)
                    title(['u(t,x,z), t=', num2str(j*dt)]);
                    xlabel('x');
                    ylabel('z');
                end
                M(j+1) = getframe(fig);
            end
        end
        
        function lambda_c = critical_lambda(initializer, pol, m, param,...
                                            lambda_name, lambda_l, lambda_r,...
                                            MaxIter, lambdaTol, betaTol)
        %critical_lambda finds the value of lambda that makes max_k Real( beta_m,k) = 0
        %   We use the Illinois algorithm. 
        %   By default lambdaTol = 1e-10, fTol = 1e-10, MaxIter = 20
            if ~exist('MaxIter', 'var')
                MaxIter = 50;
            end
            if ~exist('lambdaTol', 'var')
                lambdaTol = 1e-10;
            end
            if ~exist('betaTol', 'var')
                betaTol = 1e-10;
            end
            param = setfield(param, lambda_name, lambda_l);
            gal = initializer(param, m, pol);
            beta_l = real(gal.eigenvalues(1, 'descend'));
            param = setfield(param, lambda_name, lambda_r);
            gal = initializer(param, m, pol);
            beta_r = real(gal.eigenvalues(1, 'descend'));
            if beta_l*beta_r>0
                error('beta(lambda_l) and beta(lambda_r) must have different signs');
            elseif abs(beta_r-beta_l)<betaTol || abs(lambda_r-lambda_l)<lambdaTol
                lambda_c = (lambda_l+lambda_r)/2;
                return;
            end
                
            for n=1:MaxIter
                lambda_c = (beta_l*lambda_r-beta_r*lambda_l)/(beta_l-beta_r);
                param = setfield(param, lambda_name, lambda_c);
                gal = initializer(param, m, pol);
                beta = real(gal.eigenvalues(1, 'descend'));
                if abs(beta_r-beta_l)<betaTol || abs(lambda_r-lambda_l)<lambdaTol
                    lambda_c = (lambda_l+lambda_r)/2;
                    return;
                elseif beta*beta_r>0
                    lambda_r = lambda_c;
                    beta_r = beta;
                elseif beta_l*beta>0
                    lambda_l = lambda_c;
                    beta_l = beta;
                else
                    break;
                end
            end
                        
        end
        
        function [lambda_2_c, m_c_] = neutral_stability_curves(initializer, pol, m_, param, lambda_1_name, lambda_1_,...  
                                               lambda_2_name, lambda_2_l, lambda_2_r, ...
                                               MaxIter, lambdaTol, betaTol)
            if ~exist('MaxIter', 'var')
                MaxIter = 50;
            end
            if ~exist('lambdaTol', 'var')
                lambdaTol = 1e-10;
            end
            if ~exist('betaTol', 'var')
                betaTol = 1e-10;
            end
            lambda_2_c = zeros(numel(m_), numel(lambda_1_));
            legendInfo = cell(numel(m_), 1);
            for i=1:numel(m_)
                m = m_(i);
                legendInfo{i} = ['m = ', num2str(m)];
                for j=1:numel(lambda_1_)
                    lambda_1 = lambda_1_(j);
                    param = setfield(param, lambda_1_name, lambda_1);
                    lambda_2_c(i,j) = TransitionUtils.critical_lambda(initializer, pol, m, param,...
                                                      lambda_2_name, lambda_2_l, lambda_2_r,...
                                                      MaxIter, lambdaTol, betaTol);
                end
                plot(lambda_1_, lambda_2_c(i,:));
                hold on
            end
            
            xlabel(lambda_1_name);
            ylabel(lambda_2_name);
            legend(legendInfo);
            
            [lambda_2_c, m_c_] = min(lambda_2_c, [], 1);
            
        end
    end
end

