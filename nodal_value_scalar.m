function U = nodal_value_scalar(Mesh,U_aux,h)

nElements = size(Mesh.Elements,1);
U = zeros(nElements,1);
for i = 1:1/h/2
    for j = 1:1/h/2
        Elements = zeros(1,4);
        ref_no = 1+2*(j-1)+(2/h+2)*(i-1);
        Dof = Mesh.Vert2Edge(ref_no+1/h+2,ref_no+1/h+1);
        Elements(1:2) = Mesh.Edge2Elem(Dof,:);
        Dof = Mesh.Vert2Edge(ref_no+1/h+2,ref_no+1/h+2+1);
        Elements(3:4) = Mesh.Edge2Elem(Dof,:);
        U(Elements) = U_aux(i+(j-1)/h/2);
    end
end