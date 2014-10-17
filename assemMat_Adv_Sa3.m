function varargout = assemMat_Adv_Sa3(Mesh,U,Dofs,h,varargin)
% Advection matrix for the saturation equation (piecewise constant)

%   Chak Shing Lee
%   Department of Mathematics
%   Texas A&M University
%   2013 July

% Initialize constants

nElements = size(Mesh.Elements,1);     % Number of elements in new mesh
% Allocate memory

I = zeros(4*nElements/4,1);
J = zeros(4*nElements/4,1);
A = zeros(4*nElements/4,1);

% Assemble element contributions

loc = 1:4;
for i = 1:1/h/2
    for j = 1:1/h/2
        J_idx = zeros(4,1);
        Aloc = zeros(4,1);
        
        idx = i+(j-1)/h/2;
        ref_no = 1+2*(j-1)+(2/h+2)*(i-1);      % bottom left corner vertice of a 2 by 2 block
        
        % loop over all edges of the current (BIG) element
        
        if j == 1
            J_idx(1) = idx;
            Aloc(1) = 0;
        else
            eidx = Mesh.Vert2Edge(ref_no,ref_no+1/h+1);
            %Element = Mesh.Edge2Elem(eidx,1);
            U_dot_n = -U(Dofs(eidx));
            eidx = Mesh.Vert2Edge(ref_no+2/h+2,ref_no+1/h+1);
            %Element = Mesh.Edge2Elem(eidx,1);
            U_dot_n = U_dot_n-U(Dofs(eidx));
            if U_dot_n > 0
                J_idx(1) = idx;
            else
                J_idx(1) = idx-1/h/2;
            end
            Aloc(1) = U_dot_n*h;
        end
        
        if i == 1/h/2
            J_idx(2) = idx;
            Aloc(2) = 0;
        else
            eidx = Mesh.Vert2Edge(ref_no+2/h+2,ref_no+2/h+2+1);
            %Element = Mesh.Edge2Elem(eidx,1);
            U_dot_n = U(Dofs(eidx));
            eidx = Mesh.Vert2Edge(ref_no+2/h+2+1,ref_no+2/h+2+2);
            %Element = Mesh.Edge2Elem(eidx,1);
            U_dot_n = U_dot_n+U(Dofs(eidx));
            if U_dot_n > 0
                J_idx(2) = idx;
            else
                J_idx(2) = idx+1;
            end
            Aloc(2) = U_dot_n*h;
        end
        
        if j == 1/h/2
            J_idx(3) = idx;
            Aloc(3) = 0;
        else
            eidx = Mesh.Vert2Edge(ref_no+2/h+2+2,ref_no+1/h+1+2);
            %Element = Mesh.Edge2Elem(eidx,1);
            U_dot_n = U(Dofs(eidx));
            eidx = Mesh.Vert2Edge(ref_no+1/h+1+2,ref_no+2);
            %Element = Mesh.Edge2Elem(eidx,1);
            U_dot_n = U_dot_n+U(Dofs(eidx));
            if U_dot_n > 0
                J_idx(3) = idx;
            else
                J_idx(3) = idx+1/h/2;
            end
            Aloc(3) = U_dot_n*h;
        end
        
        if i == 1
            J_idx(4) = idx;
            Aloc(4) = 0;
        else
            eidx = Mesh.Vert2Edge(ref_no+2,ref_no+1);
            %Element = Mesh.Edge2Elem(eidx,1);
            U_dot_n = -U(Dofs(eidx));
            eidx = Mesh.Vert2Edge(ref_no,ref_no+1);
            %Element = Mesh.Edge2Elem(eidx,1);
            U_dot_n = U_dot_n-U(Dofs(eidx));
            if U_dot_n > 0
                J_idx(4) = idx;
            else
                J_idx(4) = idx-1;
            end
            Aloc(4) = U_dot_n*h;
        end
        
        % Add contributions to stiffness matrix
        
        I(loc) = set_Rows(idx,4);
        J(loc) = set_Cols(J_idx',1);
        A(loc) = Aloc(:);
        loc = loc+4;
    end
end

% Assign output arguments

if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I,J,A);
end

return