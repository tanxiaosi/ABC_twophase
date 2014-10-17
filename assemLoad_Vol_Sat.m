function L = assemLoad_Vol_Sat(Mesh,FHandle,QuadRule,varargin)
% Load vector for the saturation equation (piecewise constant)

%   Chak Shing Lee
%   Department of Mathematics
%   Texas A&M University
%   2013 July

% Initialize constants

nElements = size(Mesh.Elements,1);     % Number of elements

N = shap_BFE(QuadRule.x);
grad_N = grad_shap_BFE(QuadRule.x);  

% Preallocate memory

L = zeros(nElements,1);

% Assemble element contributions

for i = 1:nElements
    
    % Extract vertex numbers
    
    vidx = Mesh.Elements(i,:);
    
    % Extract vertices and tangents of current element
    
    P1 = Mesh.Coordinates(vidx(1),:);
    P2 = Mesh.Coordinates(vidx(2),:);
    P3 = Mesh.Coordinates(vidx(3),:);
    P4 = Mesh.Coordinates(vidx(4),:);
    
    x = N(:,1)*P1+N(:,2)*P2+N(:,3)*P3+N(:,4)*P4;
    z1 = P1(1)*grad_N(:,1:2)+P2(1)*grad_N(:,3:4) + ...
         P3(1)*grad_N(:,5:6)+P4(1)*grad_N(:,7:8);
    z2 = P1(2)*grad_N(:,1:2)+P2(2)*grad_N(:,3:4) + ...
         P3(2)*grad_N(:,5:6)+P4(2)*grad_N(:,7:8);
    det_DPhi_K = abs(z1(:,1).*z2(:,2)-z1(:,2).*z2(:,1));
    
    % Compute function values
    
    FVal = FHandle(x,varargin{:});
    
    Lloc = sum(QuadRule.w.*FVal.*det_DPhi_K);
    
    % Add contributions to stiffness matrix
    
    L(i,1) = L(i,1)+Lloc;
    
end

return