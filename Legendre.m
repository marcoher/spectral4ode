classdef Legendre
    %Legendre polynomials
    %   Detailed explanation goes here
    
    properties
        N
        tol
        x
        w
        L
        dL
        d2L
        d3L
        d4L
        norm2
    end
    
    methods
        function obj=Legendre(N,tol)
            obj.N=N;
            obj.tol=tol;
            obj.x=zeros(N+1,1);
            obj.w=zeros(N+1,1);
            obj.L=zeros(N+1);
            obj.dL=zeros(N+1);
            obj.d2L=zeros(N+1);
            obj.d3L=zeros(N+1);
            obj.d4L=zeros(N+1);
            obj.norm2=zeros(N+1,1);
            obj.x(1)=-1; 
            obj.x(N+1)=1;
            obj.w(1)=2/(N*(N+1)*dm_Leg_n(N,0,-1)^2); 
            obj.w(N+1)=2/(N*(N+1)*dm_Leg_n(N,0,1)^2);
            for k=0:N
                obj.L(k+1,1)=dm_Leg_n(k,0,-1);
                obj.dL(k+1,1)=dm_Leg_n(k,1,-1);
                obj.d2L(k+1,1)=dm_Leg_n(k,2,-1);
                obj.d3L(k+1,1)=dm_Leg_n(k,3,-1);
                obj.d4L(k+1,1)=dm_Leg_n(k,4,-1);
                obj.L(k+1,N+1)=dm_Leg_n(k,0,1);
                obj.dL(k+1,N+1)=dm_Leg_n(k,1,1);
                obj.d2L(k+1,N+1)=dm_Leg_n(k,2,1);
                obj.d3L(k+1,N+1)=dm_Leg_n(k,3,1);
                obj.d4L(k+1,N+1)=dm_Leg_n(k,4,1);
                obj.norm2(k+1)=2/(2*k+1);
            end
            for j=1:N-1
                theta1=(4*j-1)*pi/(4*N+2);
                theta2=(4*(j+1)-1)*pi/(4*N+2);
                sigma1=-(1-(N-1)/(8*N^3)...
                    -1/(384*N^4)*(39-28/sin(theta1)^2))*cos(theta1);
                sigma2=-(1-(N-1)/(8*N^3)...
                    -1/(384*N^4)*(39-28/sin(theta2)^2))*cos(theta2);
                z=(sigma1+sigma2)/2;
                error=1;
                while error>tol
                    znew=z-(1-z^2)*dm_Leg_n(N,1,z)/...
                        (2*z*dm_Leg_n(N,1,z)-N*(N+1)*dm_Leg_n(N,0,z));
                    error=abs(z-znew);
                    z=znew;
                end
                obj.x(j+1)=z;
                obj.w(j+1)=2/(N*(N+1)*dm_Leg_n(N,0,z)^2);
                for k=0:N
                    obj.L(k+1,j+1)=dm_Leg_n(k,0,z);
                    obj.dL(k+1,j+1)=dm_Leg_n(k,1,z);
                    obj.d2L(k+1,j+1)=dm_Leg_n(k,2,z);
                    obj.d3L(k+1,j+1)=dm_Leg_n(k,3,z);
                    obj.d4L(k+1,j+1)=dm_Leg_n(k,4,z);
                end
            end
        end
        
        function v = FLegT( obj, u )
            v=((0:obj.N)'+0.5).*(obj.L*((obj.w).*u));
            v(obj.N+1)=obj.N/(2*obj.N+1)*v(obj.N+1);
        end

        function [u,du,d2u,d3u,d4u] = BLegT( obj, v )
            switch nargout
                case 1
                    u=obj.L'*v;
                case 2
                    u=obj.L'*v;
                    du=obj.dL'*v;
                case 3
                    u=obj.L'*v;
                    du=obj.dL'*v;
                    d3u=obj.d3L'*v;
                case 4
                    u=obj.L'*v;
                    du=obj.dL'*v;
                    d2u=obj.d2L'*v;
                    d3u=obj.d3L'*v;
                case 5
                    u=obj.L'*v;
                    du=obj.dL'*v;
                    d2u=obj.d2L'*v;
                    d3u=obj.d3L'*v;
                    d4u=obj.d4L'*v;
            end
        end
        
    end
    
end

function r = dm_Leg_n( n, m, x )
    s=zeros(n+1,m+1);
    for j=0:m
        if j==0
            s(1,j+1)=1;
            s(2,j+1)=x;
            for k=1:n-1
               s(k+2,j+1)=((2*k+1)*x*s(k+1,j+1)-k*s(k,j+1))/(k+1);
            end
        else
            s(1,j+1)=0;
            if j==1
                s(2,j+1)=1;
            else
                s(2,j+1)=0;
            end
            for k=1:n-1
                s(k+2,j+1)=(2*k+1)*s(k+1,j)+s(k,j+1);
            end
        end
    end
    r=s(n+1,m+1);
end
