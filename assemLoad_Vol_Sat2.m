function L = assemLoad_Vol_Sat2(Mesh,FHandle,QuadRule,h,varargin)
% Load vector for the saturation equation (piecewise constant)

%   Chak Shing Lee
%   Department of Mathematics
%   Texas A&M University
%   2013 July

% Initialize constants

nElements = size(Mesh.Elements,1);     % Number of elements

% Preallocate memory

L = zeros(nElements,1);

% Assemble element contributions

for i = 1:nElements
    
    % Extract vertex numbers
    
    ref_no = Mesh.Elements(i,1);
    
    x = ones(size(QuadRule.x,1),1)*Mesh.Coordinates(ref_no,:)+QuadRule.x*h;

    
    % Compute function values
    
    FVal = FHandle(x,varargin{:});
    
    Lloc = sum(QuadRule.w.*FVal)*h^2;
    
    % Add contributions to stiffness matrix
    
    L(i,1) = L(i,1)+Lloc;
    
end

return