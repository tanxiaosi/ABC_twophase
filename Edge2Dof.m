function Dofs = Edge2Dof(Mesh,mark)
Dofs = zeros(size(Mesh.Edges,1),1);
for i = 1:size(mark,2)
    Dofs(mark(i)) = i;
end