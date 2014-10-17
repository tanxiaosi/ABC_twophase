function varargout = assemMat_Adv_Sat(Mesh,U,varargin)
% Advection matrix for the saturation equation (piecewise constant)

%   Chak Shing Lee
%   Department of Mathematics
%   Texas A&M University
%   2013 July

% Initialize constants

nElements = size(Mesh.Elements,1);     % Number of elements in new mesh


% Allocate memory

I = zeros(4*nElements,1);
J = zeros(4*nElements,1);
A = zeros(4*nElements,1);

% Assemble element contributions

loc = 1:4;
for i = 1:nElements
    
    J_idx = zeros(4,1);
    Aloc = zeros(4,1);
    
    % Extract vertices and velocity value of current element
    
    vidx = Mesh.Elements(i,:);
    vertices = Mesh.Coordinates(vidx,:);
    normals = [vertices(2,:)-vertices(1,:);vertices(3,:)-vertices(2,:);vertices(4,:)-vertices(3,:);vertices(1,:)-vertices(4,:)]*[0 -1;1 0];
    normals(1,:) = normals(1,:)/norm(normals(1,:));
    normals(2,:) = normals(2,:)/norm(normals(2,:));
    normals(3,:) = normals(3,:)/norm(normals(3,:));
    normals(4,:) = normals(4,:)/norm(normals(4,:));
    
    % loop over all edges of the current element
    
    eidx = Mesh.Vert2Edge(vidx(1),vidx(2));
    if Mesh.BdFlags(eidx) < 0
        J_idx(1) = i;
        Aloc(1) = 0;
    else
        elem_adj = setdiff(Mesh.Edge2Elem(eidx,:),i);
        U_dot_n = normals(1,:)*(U([2*i-1 2*i])+U([2*elem_adj-1 2*elem_adj]))/2;
        if U_dot_n > 0 
            J_idx(1) = i;
        else
            J_idx(1) = elem_adj;
        end
        Aloc(1) = U_dot_n*norm(vertices(2,:)-vertices(1,:));
    end
    
    eidx = Mesh.Vert2Edge(vidx(3),vidx(2));
    if Mesh.BdFlags(eidx) < 0
        J_idx(2) = i;
        Aloc(2) = 0;
    else
        elem_adj = setdiff(Mesh.Edge2Elem(eidx,:),i);
        U_dot_n = normals(2,:)*(U([2*i-1 2*i])+U([2*elem_adj-1 2*elem_adj]))/2;
        if U_dot_n > 0 
            J_idx(2) = i;
        else
            J_idx(2) = elem_adj;
        end
        Aloc(2) = U_dot_n*norm(vertices(2,:)-vertices(3,:));
    end
    
    eidx = Mesh.Vert2Edge(vidx(3),vidx(4));
    if Mesh.BdFlags(eidx) < 0
        J_idx(3) = i;
        Aloc(3) = 0;
    else
        elem_adj = setdiff(Mesh.Edge2Elem(eidx,:),i);
        U_dot_n = normals(3,:)*(U([2*i-1 2*i])+U([2*elem_adj-1 2*elem_adj]))/2;
        if U_dot_n > 0 
            J_idx(3) = i;
        else
            J_idx(3) = elem_adj;
        end
        Aloc(3) = U_dot_n*norm(vertices(4,:)-vertices(3,:));
    end
    
    eidx = Mesh.Vert2Edge(vidx(1),vidx(4));
    if Mesh.BdFlags(eidx) < 0
        J_idx(4) = i;
        Aloc(4) = 0;
    else
        elem_adj = setdiff(Mesh.Edge2Elem(eidx,:),i);
        U_dot_n = normals(4,:)*(U([2*i-1 2*i])+U([2*elem_adj-1 2*elem_adj]))/2;
        if U_dot_n > 0 
            J_idx(4) = i;
        else
            J_idx(4) = elem_adj;
        end
        Aloc(4) = U_dot_n*norm(vertices(4,:)-vertices(1,:));
    end
    
    % Add contributions to stiffness matrix
    
    I(loc) = set_Rows(i,4);
    J(loc) = set_Cols(J_idx',1);
    A(loc) = Aloc(:);
    loc = loc+4;
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