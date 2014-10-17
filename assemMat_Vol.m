function M_u = assemMat_Vol(Mesh,h,mark,coef,S,mu_w,mu_o)
n=1/h;
m = zeros(n^2-2*n,1);
I = 1:(n^2-2*n);
J = 1:(n^2-2*n);
for i = 1:size(mark,2)
    Elements = Mesh.Edge2Elem(mark(i),:);
    if Elements(1)~=0
        lambda_S = S(Elements(1))^2/mu_w + (S(Elements(1))-1)^2/mu_o;
        m(i)=h^2/(lambda_S*coef(Elements(1)));
    end
    if Elements(2)~=0
        lambda_S = S(Elements(2))^2/mu_w + (S(Elements(2))-1)^2/mu_o;
        m(i) = m(i) + h^2/(lambda_S*coef(Elements(2)));
    end
end

M_u =sparse(I,J,m);