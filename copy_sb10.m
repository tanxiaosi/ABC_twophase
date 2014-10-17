function [kappa,kappa2] = copy_sb10(Mesh,N)

load Udata.mat
K0=1e14*KU(:,:,:,85);

% Initialize constants

nElements = size(Mesh.Elements,1);     % Number of elements

% Generate the value of permeability at each vertex nu

kappa = zeros(nElements,1);
kappa2 = zeros(nElements,1);

for i = 1:N*(N-1)*2
    Edges = find(Mesh.DDMData(1,:)==i);
    for j = 1:numel(Edges)
        Elem = Mesh.Edge2Elem(Edges(j),:);
        if Elem(1) ~= 0 && Elem(2) ~= 0
            kappa2(Elem(1)) = 1;
            kappa2(Elem(2)) = 1;
        end
    end
end
for i =1:nElements
    vertices = sum(Mesh.Coordinates(Mesh.Elements(i,:),:))/4;
    x = ceil(vertices(1)*60);
    x = x + (x==0);
    y = ceil(vertices(2)*220);
    y = y + (y==0);
    kappa(i) = K0(1,x,y);
    kappa2(i) = kappa2(i)*kappa(i);
end