function kappa = perm_k(K_Handle,Mesh)

% Initialize constants

nElements = size(Mesh.Elements,1);     % Number of elements

% Generate the value of permeability in each fine grid element

kappa = zeros(nElements,1);

for i =1:nElements
    kappa(i) = K_Handle(sum(Mesh.Coordinates(Mesh.Elements(i,:),:))/4);
end