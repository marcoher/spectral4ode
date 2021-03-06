classdef Chebyshev
    %Chebyshev Class for Chebyshev polynomials
    %   Contains basic methods for constructing Chebyshev polynomials using
    %   Chebyshev-Gauss-Lobatto quadrature, finding nodes, weigths, and
    %   forward/backward transforms.
    
    properties
        N       % Max degree [10]
        x       % LGL Nodes x_0,...,x_N: (N+1)x1
        w       % Weigths w_0,...,w_N: (N+1)x1
        P       % Polynomial values P_i(x_j): (N+1)x(N+1)
        dP      % Polynomial values P_i'(x_j): (N+1)x(N+1)
        d2P     % Polynomial values P_i''(x_j): (N+1)x(N+1)
        d3P     % Polynomial values P_i'''(x_j): (N+1)x(N+1)
        d4P     % Polynomial values P_i''''(x_j): (N+1)x(N+1)
        norm2   % L^2 norm of P_i: (N+1)x1
    end
    
    methods
        function obj=Chebyshev(N)
            if nargin==0
                N=10;
            end
            obj.N=N;
                     
            obj.x = cos(pi*(N:-1:0)'/N);
            
            obj.w = ones(N+1,1)*pi/2;
            obj.w(1) = pi/4; 
            obj.w(N+1) = pi/4;
            
            obj.P = zeros(N+1);
            obj.dP = zeros(N+1);
            obj.d2P = zeros(N+1);
            obj.d3P = zeros(N+1);
            obj.d4P = zeros(N+1);
            
            obj.norm2 = ones(N+1,1)*pi/2;
            obj.norm2(1) = pi;
                        
            for k=1:numel(obj.x)
                [obj.P(:,k), obj.dP(:,k), obj.d2P(:,k),...
                    obj.d3P(:,k), obj.d4P(:,k)] = Cheb(N,obj.x(k));
            end
            
        end
        
        function v = FT( obj, u )
            % Forward transform
            vv = ifft([u(obj.N+1:-1:1);u(2:obj.N)]);
            v = [vv(1); 2*vv(2:obj.N); vv(obj.N+1)];
        end

        function [u,du,d2u,d3u,d4u] = BT( obj, v )
            % Backward transform ... should speed up with fft?
            switch nargout
                case 1
                    u=obj.P'*v;
                case 2
                    u=obj.P'*v;
                    du=obj.dP'*v;
                case 3
                    u=obj.P'*v;
                    du=obj.dP'*v;
                    d3u=obj.d3P'*v;
                case 4
                    u=obj.P'*v;
                    du=obj.dP'*v;
                    d2u=obj.d2P'*v;
                    d3u=obj.d3P'*v;
                case 5
                    u=obj.P'*v;
                    du=obj.dP'*v;
                    d2u=obj.d2P'*v;
                    d3u=obj.d3P'*v;
                    d4u=obj.d4P'*v;
            end
        end
        
    end
    
end

function [T, dT, d2T, d3T, d4T] = Cheb( n, x )
    m=4;
    s=zeros(n+1,m+1);
    for j=0:m
        if j==0
            s(1,j+1)=1;
            s(2,j+1)=x;
            for k=1:n-1
               s(k+2,j+1)=2*x*s(k+1,j+1)-s(k,j+1); 
            end
        else
            s(1,j+1)=0;
            if j==1
                s(2,j+1)=1;
            else
                s(2,j+1)=0;
            end
            for k=1:n-1
                s(k+2,j+1) = 2*x*s(k+1,j+1)+2*j*s(k+1,j)-s(k,j+1);
            end
        end
    end
    T=s(:,m+1-4);
    dT=s(:,m+1-3);
    d2T=s(:,m+1-2);
    d3T=s(:,m+1-1);
    d4T=s(:,m+1);
    
end
