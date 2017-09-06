classdef Galerkin
    %Galerkin Class for spectral Galerkin representation of ODE
    %   Contains basic methods for converting between frequency/physical
    %   space, finding eigenfunctions/eigenvalues, plotting the result
    %   assuming axial symmetry, and solving linear systems. The underlying
    %   basis of orthogonal polynomials must be set a priori.
    
    properties
        pol     % Basis of orthogonal polynomials (Legendre or Chebyshev)
        neq     % Number of equations
        ord     % Order of each equation (either 4 or 2)
        coeff   % Coefficients in the Galerkin decomposition
        opRhs   % Differential operator of order 4 or 2
        opLhs   % Differential operator of order 2 or 0
        x_fun   % Coordinate transform (linear, increasing)
        sc      % Scale factor
    end
    
    methods
        function obj = Galerkin(pol, bc_l, bc_r, rhs, lhs, x_l, x_r)
            
            obj.neq = numel(bc_l);
            
            if ~(numel(bc_r)==obj.neq)
                error(['Number of boundary conditions on the left does '...
                    'not match number of boundary conditions on the right']);
            end
            
            obj.pol = pol;
            
            obj.sc = 2.0/(x_r - x_l);
            
            obj.x_fun = @(x) (1 + x)/2 * x_r + (1 - x)/2 * x_l;
            
            [obj.ord, obj.coeff] = obj.basis_expansion(bc_l, bc_r);

            [obj.opLhs, obj.opRhs] = obj.construct_ops(lhs, rhs);
            
        end
            
        
        function [ord, coeff] = basis_expansion(obj, bc_l, bc_r)
            coeff = cell(obj.neq,1);
            ord = zeros(obj.neq,1);
            
            for i=1:obj.neq
                ord(i) = 2*numel(bc_l{i});
                if ~(2*numel(bc_r{i})==ord(i))
                    error(['%d variable has different number of boundary'...
                        ' conditions one the left and on the right'],i);
                end
                switch ord(i)
                    case 4
                        A = zeros(4); 
                        A0 = zeros(4,1);
                        c_ = zeros(obj.pol.N-3,4);
                        for p=0:obj.pol.N-4
                            for j=1:4
                                A(1,j) = bc_l{i}{1}(...
                                    obj.pol.P(p+1+j,1),...
                                    obj.sc*obj.pol.dP(p+1+j,1),...
                                    obj.sc^2*obj.pol.d2P(p+1+j,1),...
                                    obj.sc^3*obj.pol.d3P(p+1+j,1) );
                                
                                A(2,j) = bc_l{i}{2}(...
                                    obj.pol.P(p+1+j,1),...
                                    obj.sc*obj.pol.dP(p+1+j,1),...
                                    obj.sc^2*obj.pol.d2P(p+1+j,1),...
                                    obj.sc^3*obj.pol.d3P(p+1+j,1) );
                                
                                A(3,j) = bc_r{i}{1}(...
                                    obj.pol.P(p+1+j,end),...
                                    obj.sc*obj.pol.dP(p+1+j,end),...
                                    obj.sc^2*obj.pol.d2P(p+1+j,end),...
                                    obj.sc^3*obj.pol.d3P(p+1+j,end) );
                                
                                A(4,j) = bc_r{i}{2}(...
                                    obj.pol.P(p+1+j,end),...
                                    obj.sc*obj.pol.dP(p+1+j,end),...
                                    obj.sc^2*obj.pol.d2P(p+1+j,end),...
                                    obj.sc^3*obj.pol.d3P(p+1+j,end) );
                            end
                            
                            A0(1) = bc_l{i}{1}(...
                                obj.pol.P(p+1,1),...
                                obj.sc*obj.pol.dP(p+1,1),...
                                obj.sc^2*obj.pol.d2P(p+1,1),...
                                obj.sc^3*obj.pol.d3P(p+1,1) );
                            
                            A0(2) = bc_l{i}{2}(...
                                obj.pol.P(p+1,1),...
                                obj.sc*obj.pol.dP(p+1,1),...
                                obj.sc^2*obj.pol.d2P(p+1,1),...
                                obj.sc^3*obj.pol.d3P(p+1,1) );
                            
                            A0(3) = bc_r{i}{1}(...
                                obj.pol.P(p+1,end),...
                                obj.sc*obj.pol.dP(p+1,end),...
                                obj.sc^2*obj.pol.d2P(p+1,end),...
                                obj.sc^3*obj.pol.d3P(p+1,end) );
                            
                            A0(4) = bc_r{i}{2}(...
                                obj.pol.P(p+1,end),...
                                obj.sc*obj.pol.dP(p+1,end),...
                                obj.sc^2*obj.pol.d2P(p+1,end),...
                                obj.sc^3*obj.pol.d3P(p+1,end) );
                            
                            c_(p+1,:)=-(A\A0)';
                        end
                    case 2
                        A = zeros(2); 
                        A0 = zeros(2,1);
                        c_ = zeros(obj.pol.N-1,2);
                        for p=0:obj.pol.N-2
                            for j=1:2
                                A(1,j) = bc_l{i}{1}(...
                                    obj.pol.P(p+1+j,1),...
                                    obj.sc*obj.pol.dP(p+1+j,1));
                                
                                A(2,j) = bc_r{i}{1}(...
                                    obj.pol.P(p+1+j,end),...
                                    obj.sc*obj.pol.dP(p+1+j,end));
                            end
                            
                            A0(1) = bc_l{i}{1}(...
                                obj.pol.P(p+1,1),...
                                obj.sc*obj.pol.dP(p+1,1));
                            
                            A0(2)=bc_r{i}{1}(...
                                obj.pol.P(p+1,end),...
                                obj.sc*obj.pol.dP(p+1,end));
                            
                            c_(p+1,:)=-(A\A0)';
                        end
                    otherwise
                        error('Equation %d must be order 2 or 4',i);
                end
                coeff{i}=c_;
            end
            
        end
        
        function [opLhs, opRhs] = construct_ops(obj, lhs, rhs)
                        
            op_rhs = cell(obj.neq, obj.neq);
            
            op_lhs = cell(obj.neq, obj.neq);
            
            for i=1:obj.neq
                for j=1:obj.neq
                    L_ = zeros(obj.pol.N-obj.ord(i)+1,...
                               obj.pol.N-obj.ord(j)+1);
                    
                    A_ = zeros(obj.pol.N-obj.ord(i)+1,...
                               obj.pol.N-obj.ord(j)+1);
                    
                    for q=0:obj.pol.N-obj.ord(j)
                        phi_qj=( obj.pol.P(q+1,:)...
                                +obj.coeff{j}(q+1,1)*obj.pol.P(q+2,:)...
                                +obj.coeff{j}(q+1,2)*obj.pol.P(q+3,:) )';
                        
                        dphi_qj=( obj.pol.dP(q+1,:)...
                                +obj.coeff{j}(q+1,1)*obj.pol.dP(q+2,:)...
                                +obj.coeff{j}(q+1,2)*obj.pol.dP(q+3,:) )';
                        
                        d2phi_qj=( obj.pol.d2P(q+1,:)...
                                +obj.coeff{j}(q+1,1)*obj.pol.d2P(q+2,:)...
                                +obj.coeff{j}(q+1,2)*obj.pol.d2P(q+3,:) )';
                        
                        if obj.ord(j)==4
                            phi_qj=phi_qj+...
                                (obj.coeff{j}(q+1,3)*obj.pol.P(q+4,:)...
                                +obj.coeff{j}(q+1,4)*obj.pol.P(q+5,:) )';
                            dphi_qj=dphi_qj+...
                                (obj.coeff{j}(q+1,3)*obj.pol.dP(q+4,:)...
                                +obj.coeff{j}(q+1,4)*obj.pol.dP(q+5,:) )';
                            d2phi_qj=d2phi_qj+...
                                (obj.coeff{j}(q+1,3)*obj.pol.d2P(q+4,:)...
                                +obj.coeff{j}(q+1,4)*obj.pol.d2P(q+5,:) )';
                            d3phi_qj=( obj.pol.d3P(q+1,:)...
                                +obj.coeff{j}(q+1,1)*obj.pol.d3P(q+2,:)...
                                +obj.coeff{j}(q+1,2)*obj.pol.d3P(q+3,:)...
                                +obj.coeff{j}(q+1,3)*obj.pol.d3P(q+4,:)...
                                +obj.coeff{j}(q+1,4)*obj.pol.d3P(q+5,:) )';
                            d4phi_qj=( obj.pol.d4P(q+1,:)...
                                +obj.coeff{j}(q+1,1)*obj.pol.d4P(q+2,:)...
                                +obj.coeff{j}(q+1,2)*obj.pol.d4P(q+3,:)...
                                +obj.coeff{j}(q+1,3)*obj.pol.d4P(q+4,:)...
                                +obj.coeff{j}(q+1,4)*obj.pol.d4P(q+5,:) )';
                            rhs_ij_phi_qj=rhs{i}{j}(...
                                obj.x_fun(obj.pol.x),...
                                phi_qj,...
                                obj.sc*dphi_qj,...
                                obj.sc^2*d2phi_qj,...
                                obj.sc^3*d3phi_qj,...
                                obj.sc^4*d4phi_qj);
                            lhs_ij_phi_qj=lhs{i}{j}(...
                                obj.x_fun(obj.pol.x),...
                                phi_qj,...
                                obj.sc*dphi_qj,...
                                obj.sc^2*d2phi_qj);
                        elseif obj.ord(j)==2
                             rhs_ij_phi_qj=rhs{i}{j}(...
                                obj.x_fun(obj.pol.x),...
                                phi_qj,...
                                obj.sc*dphi_qj,...
                                obj.sc^2*d2phi_qj);
                            lhs_ij_phi_qj=lhs{i}{j}(...
                                obj.x_fun(obj.pol.x),...
                                phi_qj);
                        end
                        rr=obj.pol.FT(rhs_ij_phi_qj);
                        ll=obj.pol.FT(lhs_ij_phi_qj);
                        if obj.ord(i)==4
                            L_(:,q+1)=obj.pol.norm2(1:end-4).*rr(1:end-4)...
                                +obj.coeff{i}(:,1).*obj.pol.norm2(2:end-3).*rr(2:end-3)...
                                +obj.coeff{i}(:,2).*obj.pol.norm2(3:end-2).*rr(3:end-2)...
                                +obj.coeff{i}(:,3).*obj.pol.norm2(4:end-1).*rr(4:end-1)...
                                +obj.coeff{i}(:,4).*obj.pol.norm2(5:end).*rr(5:end);
                            A_(:,q+1)=obj.pol.norm2(1:end-4).*ll(1:end-4)...
                                +obj.coeff{i}(:,1).*obj.pol.norm2(2:end-3).*ll(2:end-3)...
                                +obj.coeff{i}(:,2).*obj.pol.norm2(3:end-2).*ll(3:end-2)...
                                +obj.coeff{i}(:,3).*obj.pol.norm2(4:end-1).*ll(4:end-1)...
                                +obj.coeff{i}(:,4).*obj.pol.norm2(5:end).*ll(5:end);
                        elseif obj.ord(i)==2
                            L_(:,q+1)=obj.pol.norm2(1:end-2).*rr(1:end-2)...
                                +obj.coeff{i}(:,1).*obj.pol.norm2(2:end-1).*rr(2:end-1)...
                                +obj.coeff{i}(:,2).*obj.pol.norm2(3:end).*rr(3:end);
                            A_(:,q+1)=obj.pol.norm2(1:end-2).*ll(1:end-2)...
                                +obj.coeff{i}(:,1).*obj.pol.norm2(2:end-1).*ll(2:end-1)...
                                +obj.coeff{i}(:,2).*obj.pol.norm2(3:end).*ll(3:end);
                        end
                    end
                    op_lhs{i,j}=A_;
                    op_rhs{i,j}=L_;
                end
            end
            opRhs=cell2mat(op_rhs);
            opLhs=cell2mat(op_lhs);
            
        end
        
        function v = FPhiT( obj, u )
        %FPhiT Foward transform    
            if ~(numel(u)==obj.neq)
                error(['Number of components does not match '...
                    'Galerkin decomposition']);
            end
            v=cell(obj.neq,1);
            for i=1:obj.neq
               vv=obj.pol.FT(u{i});
               switch size(obj.coeff{i},2)
                   case 4
                       v{i}=zeros(obj.pol.N-3,1);
                       v{i}(1)=vv(1);
                       v{i}(2)=vv(2)-obj.coeff{i}(1,1)*vv(1);
                       v{i}(3)=vv(3)-obj.coeff{i}(2,1)*vv(2) ...
                           -obj.coeff{i}(1,2)*vv(1);
                       v{i}(4)=vv(4)-obj.coeff{i}(3,1)*vv(3) ...
                           -obj.coeff{i}(2,2)*vv(2) ...
                           -obj.coeff{i}(1,3)*vv(1);
                       for k=4:obj.pol.N-4
                          v{i}(k+1)=vv(k+1) ...
                              -obj.coeff{i}(k,1)*v{i}(k) ...
                              -obj.coeff{i}(k-1,2)*v{i}(k-1) ...
                              -obj.coeff{i}(k-2,3)*v{i}(k-2) ...
                              -obj.coeff{i}(k-3,4)*v{i}(k-3) ;
                       end
                   case 2
                       v{i}=zeros(obj.pol.N-1,1);
                       v{i}(1)=vv(1);
                       v{i}(2)=vv(2)-obj.coeff{i}(1,1)*vv(1);
                       for k=2:obj.pol.N-2
                          v{i}(k+1)=vv(k+1)-obj.coeff{i}(k,1)*v{i}(k)-obj.coeff{i}(k-1,2)*v{i}(k-1); 
                       end
                   otherwise
                       error('Unexpected size of c(%d,:)',i);
               end
            end
            v=cell2mat(v);
        end
        
        function [u,du,d2u,d3u] = BPhiT( obj, v )
        %BPhiT Backward transform   
            Ni=@(i) (obj.pol.N-1)*(size(obj.coeff{i},2)==2)...
                        +(obj.pol.N-3)*(size(obj.coeff{i},2)==4);
            switch nargout
                case 1
                    u=cell(obj.neq,1);
                    n=0;
                    for i=1:obj.neq
                        if size(obj.coeff{i},2)==4
                            phi=obj.pol.P(1:end-4,:)...
                                +diag(obj.coeff{i}(:,1))*obj.pol.P(2:end-3,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.P(3:end-2,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.P(4:end-1,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.P(5:end,:);
                        elseif size(obj.coeff{i},2)==2
                            phi=obj.pol.P(1:end-2,:)...
                                +diag(obj.coeff{i}(:,1))*obj.pol.P(2:end-1,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.P(3:end,:);
                        end
                        u{i}=(v(n+1:n+Ni(i)).'*phi).';
                        n=n+Ni(i);
                    end
                case 2
                    u=cell(obj.neq,1);
                    du=cell(obj.neq,1);
                    n=0;
                    for i=1:obj.neq
                        if size(obj.coeff{i},2)==4
                            phi=obj.pol.P(1:end-4,:)...
                                +diag(obj.coeff{i}(:,1))*obj.pol.P(2:end-3,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.P(3:end-2,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.P(4:end-1,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.P(5:end,:);
                            dphi=obj.pol.dP(1:end-4,:)...
                                +diag(obj.coeff{i}(:,1))*obj.pol.dP(2:end-3,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.dP(3:end-2,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.dP(4:end-1,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.dP(5:end,:);
                        elseif size(obj.coeff{i},2)==2
                            phi=obj.pol.P(1:end-2,:)...
                                +diag(obj.coeff{i}(:,1))*obj.pol.P(2:end-1,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.P(3:end,:);
                            dphi=obj.pol.dP(1:end-2,:)...
                                +diag(obj.coeff{i}(:,1))*obj.pol.dP(2:end-1,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.dP(3:end,:);
                        end
                        u{i}=( v(n+1:n+Ni(i)).'*phi ).';
                        du{i}=obj.sc*( v(n+1:n+Ni(i)).'*dphi ).';
                        n=n+Ni(i);
                    end
                case 3
                    u=cell(obj.neq,1);
                    du=cell(obj.neq,1);
                    d2u=cell(obj.neq,1);
                    n=0;
                    for i=1:obj.neq
                        if size(obj.coeff{i},2)==4
                            phi=obj.pol.P(1:end-4,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.P(2:obj.pol.N-2,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.P(3:obj.pol.N-1,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.P(4:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.P(5:obj.pol.N+1,:);
                            dphi=obj.pol.dP(1:end-4,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.dP(2:obj.pol.N-2,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.dP(3:obj.pol.N-1,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.dP(4:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.dP(5:obj.pol.N+1,:);
                            d2phi=obj.pol.d2P(1:end-4,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.d2P(2:obj.pol.N-2,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.d2P(3:obj.pol.N-1,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.d2P(4:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.d2P(5:obj.pol.N+1,:);
                        elseif size(obj.coeff{i},2)==2
                            phi=obj.pol.P(1:obj.pol.N-1,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.P(2:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.P(3:obj.pol.N+1,:);
                            dphi=obj.pol.dP(1:obj.pol.N-1,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.dP(2:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.dP(3:obj.pol.N+1,:);
                            d2phi=obj.pol.d2P(1:obj.pol.N-1,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.d2P(2:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.d2P(3:obj.pol.N+1,:);
                        end
                        u{i}=( v(n+1:n+Ni(i)).'*phi ).';
                        du{i}=obj.sc*( v(n+1:n+Ni(i)).'*dphi ).';
                        d2u{i}=obj.sc^2*( v(n+1:n+Ni(i)).'*d2phi ).';
                        n=n+Ni(i);
                    end
                case 4
                    u=cell(obj.neq,1);
                    du=cell(obj.neq,1);
                    d2u=cell(obj.neq,1);
                    d3u=cell(obj.neq,1);
                    Ni=@(i) (obj.pol.N-1)*(size(obj.coeff{i},2)==2)...
                        +(obj.pol.N-3)*(size(obj.coeff{i},2)==4);
                    n=0;
                    for i=1:obj.neq
                        if size(obj.coeff{i},2)==4
                            phi=obj.pol.P(1:obj.pol.N-3,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.P(2:obj.pol.N-2,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.P(3:obj.pol.N-1,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.P(4:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.P(5:obj.pol.N+1,:);
                            dphi=obj.pol.dP(1:obj.pol.N-3,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.dP(2:obj.pol.N-2,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.dP(3:obj.pol.N-1,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.dP(4:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.dP(5:obj.pol.N+1,:);
                            d2phi=obj.pol.d2P(1:obj.pol.N-3,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.d2P(2:obj.pol.N-2,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.d2P(3:obj.pol.N-1,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.d2P(4:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.d2P(5:obj.pol.N+1,:);
                            d3phi=obj.pol.d3P(1:obj.pol.N-3,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.d3P(2:obj.pol.N-2,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.d3P(3:obj.pol.N-1,:)...
                                +diag(obj.coeff{i}(:,3))*obj.pol.d3P(4:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,4))*obj.pol.d3P(5:obj.pol.N+1,:);
                        elseif size(obj.coeff{i},2)==2
                            phi=obj.pol.P(1:obj.pol.N-1,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.P(2:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.P(3:obj.pol.N+1,:);
                            dphi=obj.pol.dP(1:obj.pol.N-1,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.dP(2:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.dP(3:obj.pol.N+1,:);
                            d2phi=obj.pol.d2P(1:obj.pol.N-1,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.d2P(2:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.d2P(3:obj.pol.N+1,:);
                            d3phi=obj.pol.d3P(1:obj.pol.N-1,:)+...
                                diag(obj.coeff{i}(:,1))*obj.pol.d3P(2:obj.pol.N,:)...
                                +diag(obj.coeff{i}(:,2))*obj.pol.d3P(3:obj.pol.N+1,:);
                        end
                        u{i}=( v(n+1:n+Ni(i)).'*phi ).';
                        du{i}=obj.sc*( v(n+1:n+Ni(i)).'*dphi ).';
                        d2u{i}=obj.sc^2*( v(n+1:n+Ni(i)).'*d2phi ).';
                        d3u{i}=obj.sc^3*( v(n+1:n+Ni(i)).'*d3phi ).';
                        n=n+Ni(i);
                    end
            end
            
        end
        
        function beta = eigenvalues( obj, nev, flag )
        %eigenvalues finds first nev eigenvalues    
           beta_=eig(obj.opRhs,obj.opLhs);
           [~,ind]=sort(real(beta_),flag); 
           beta=beta_(ind(1:nev));
        end
        
        function [beta,u,du,d2u,d3u] = eigenfunctions( obj, nev, flag )
        %eigenfunctions finds first nev eigenfunctions and derivatives
            if strcmp(flag,'normal') || strcmp(flag,'raw')
                [v,beta]=eig(obj.opRhs,obj.opLhs);
            elseif strcmp(flag,'dual')
                [v,beta]=eig(obj.opRhs',obj.opLhs');
            else
                error('Last argument must be ''normal'', ''dual'' or ''raw''.');
            end
            [re_beta,ind]=sort(real(diag(beta)),'descend'); 
            sorted_v=zeros(size(v)); 
            sorted_beta=zeros(size(beta));
            for j=1:length(re_beta)
                sorted_v(:,j)=v(:,ind(j));  
                sorted_beta(j,j)=beta(ind(j),ind(j));
            end
            v=sorted_v(:,1:nev); 
            beta=diag(sorted_beta(1:nev,1:nev));
            if strcmp(flag,'raw') || strcmp(flag,'dual')
                u=v;
                if nargout>2
                    error('Only one output is given in frequency form.');
                end
                return;
            end
            u=cell(nev,1);
            du=cell(nev,1);
            d2u=cell(nev,1);
            d3u=cell(nev,1);
            
            for i=1:nev
                [u{i},du{i},d2u{i},d3u{i}]=obj.BPhiT(v(:,i));
            end
        end
        
        function h = plot_as_flat( obj, u, p_vals, m, flag )
        %plot_as_flat plots Re ( sum_j u_j(x) e^{i m_j p} )
        % assumes x is vertical, p is horizontal, p in p_vals interval
           
           [pp, zz] = meshgrid(p_vals, obj.x_fun(obj.pol.x));
           uu = zeros(size(pp));
           per = p_vals(end);
           if strcmp(flag,'real')
               for j=1:length(p_vals)
                   uu(:,j) = real(u)*cos(2*pi*m/per*p_vals(j))-imag(u)*sin(2*pi*m/per*p_vals(j));
               end
           elseif strcmp(flag,'imag')
               for j=1:length(p_vals)
                   uu(:,j) = real(u)*sin(2*pi*m/per*p_vals(j))+imag(u)*cos(2*pi*m/per*p_vals(j));
               end
           else
               error('Last argument must be ''real'' or ''imag''.');
           end
           h = contourf(pp,zz,uu);
        end
        
        function h = surf_as_flat( obj, u, p_vals, m, flag )
        %plot_as_flat plots Re ( u(x) e^{i m t} )
        % assumes x is vertical, t is horizontal, t in p_vals interval
           [pp,zz]=meshgrid(p_vals,obj.x_fun(obj.pol.x));
           uu=zeros(size(pp));
           if strcmp(flag,'real')
               for j=1:length(p_vals)
                   uu(:,j)=real(u)*cos(m*p_vals(j))-imag(u)*sin(m*p_vals(j));
               end
           elseif strcmp(flag,'imag')
               for j=1:length(p_vals)
                   uu(:,j)=real(u)*sin(m*p_vals(j))+imag(u)*cos(m*p_vals(j));
               end
           else
               error('Last argument must be real or imag.');
           end
           h=surf(pp,zz,uu);
        end
        
        function h = plot_as_polar( obj, u, p_vals, m, flag )
        %plot_as_polar plots Re ( u(x) e^{i m t} )
        % assumes x is radial, t is circular (so pvals is 0:dt:2*pi)
           [pp,zz]=meshgrid(p_vals,obj.x_fun(obj.pol.x));
           uu=zeros(size(pp));
           if strcmp(flag,'real')
               for j=1:length(p_vals)
                   uu(:,j)=real(u)*cos(m*p_vals(j))-imag(u)*sin(m*p_vals(j));
               end
           elseif strcmp(flag,'imag')
               for j=1:length(p_vals)
                   uu(:,j)=real(u)*sin(m*p_vals(j))+imag(u)*cos(m*p_vals(j));
               end
           else
               error('Last argument must be ''real'' or ''imag''.');
           end
           [xx,yy,uu]=pol2cart(pp,zz,uu);
           h=contourf(xx,yy,uu);
        end
        
        function [u,du,d2u,d3u] = solve_linear(obj, k, f)
        %solve_linear solves the system kAu = Lu + f
            f_=cell(obj.neq,1);
            for i=1:obj.neq
                rr=obj.pol.FT(f{i});
                switch obj.ord(i)
                    case 4
                        f_{i}=obj.pol.norm2(1:end-4).*rr(1:end-4)...
                                +obj.coeff{i}(:,1).*obj.pol.norm2(2:end-3).*rr(2:end-3)...
                                +obj.coeff{i}(:,2).*obj.pol.norm2(3:end-2).*rr(3:end-2)...
                                +obj.coeff{i}(:,3).*obj.pol.norm2(4:end-1).*rr(4:end-1)...
                                +obj.coeff{i}(:,4).*obj.pol.norm2(5:end).*rr(5:end);
                    case 2
                        f_{i}=obj.pol.norm2(1:end-2).*rr(1:end-2)...
                                +obj.coeff{i}(:,1).*obj.pol.norm2(2:end-1).*rr(2:end-1)...
                                +obj.coeff{i}(:,2).*obj.pol.norm2(3:end).*rr(3:end);
                end
            end
            f_=cell2mat(f_);
            u_=(k*obj.opLhs-obj.opRhs)\f_;
            [u,du,d2u,d3u]=obj.BPhiT(u_);
        end
        
        function I = inner(obj, u, v)
           I=cellfun(@(z1,z2) sum(obj.pol.w.*z1.*conj(z2))/obj.sc,u,v,'UniformOutput',false);
           I=cell2mat(I);
           I=sum(I(:));
        end
        
        function I = normL2(obj, u)
            I=0;
            for i=1:obj.neq
               uu=obj.pol.FT(u{i});
               I=I+sum(abs(uu).^2.*obj.pol.norm2);
            end
            I=sqrt(I);
            
        end
    end
    
end

