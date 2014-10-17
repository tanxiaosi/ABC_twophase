function M = assemMat_Vol_scalar(Mesh,h,coef)
n=1/h;

M = zeros(n^2/4,1);
I = zeros(n^2/4,1);
J = zeros(n^2/4,1);
cnt = 1;
for i = 1:1/h/2
    for j = 1:1/h/2
        Elements = zeros(1,4);
        ref_no = 1+2*(j-1)+(2/h+2)*(i-1);
        Dof = Mesh.Vert2Edge(ref_no+1/h+2,ref_no+1/h+1);
        Elements(1:2) = Mesh.Edge2Elem(Dof,:);
        Dof = Mesh.Vert2Edge(ref_no+1/h+2,ref_no+1/h+2+1);
        Elements(3:4) = Mesh.Edge2Elem(Dof,:);
        I(cnt) = i+(j-1)/h/2;
        J(cnt) = i+(j-1)/h/2;
        M(cnt) = h^2*coef(Elements(1)) + h^2*coef(Elements(2))...
            + h^2*coef(Elements(3)) + h^2*coef(Elements(4));
        cnt = cnt + 1;
    end
end

M = sparse(I,J,M);