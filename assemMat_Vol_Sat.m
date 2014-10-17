function varargout = assemMat_Vol_Sat(Mesh,QuadRule,varargin)
% Mass matrix for the saturation equation (piecewise constant)

%   Chak Shing Lee
%   Department of Mathematics
%   Texas A&M University
%   2013 July

  % Initialize constants

  nElements = size(Mesh.Elements,1);     % Number of elements in new mesh

  
  % Allocate memory
  
  I = zeros(nElements,1);
  J = zeros(nElements,1);
  A = zeros(nElements,1);
  
  % Assemble element contributions
  
  grad_N = grad_shap_BFE(QuadRule.x);  
  loc = 1;
  for i = 1:nElements
    
    % Extract vertices and tangents of current element
    
    vidx = Mesh.Elements(i,:);
    
    % Compute element mapping  
      
    P1 = Mesh.Coordinates(vidx(1),:);
    P2 = Mesh.Coordinates(vidx(2),:);
    P3 = Mesh.Coordinates(vidx(3),:);
    P4 = Mesh.Coordinates(vidx(4),:);
    
    z1 = P1(1)*grad_N(:,1:2)+P2(1)*grad_N(:,3:4) + ...
         P3(1)*grad_N(:,5:6)+P4(1)*grad_N(:,7:8);
    z2 = P1(2)*grad_N(:,1:2)+P2(2)*grad_N(:,3:4) + ...
         P3(2)*grad_N(:,5:6)+P4(2)*grad_N(:,7:8);
    det_DPhi_K = sum(QuadRule.w.*abs(z1(:,1).*z2(:,2)-z1(:,2).*z2(:,1)));

    % Compute local integral
    
    Aloc = 1/det_DPhi_K;

    % Add contributions to stiffness matrix
    
    I(loc) = set_Rows(i,1);
    J(loc) = set_Cols(i,1);
    A(loc) = Aloc(:);
    loc = loc+1;
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