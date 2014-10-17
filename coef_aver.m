function perm_harm_aver = coef_aver(Mesh,N,h,coef,S,mu_w,mu_o)

n = 1/h;
perm_harm_aver = zeros(N*(N-1)*2/(N/n),1);
lambda_S = S.^2/mu_w + (S-1).^2/mu_o;

for kk = 1:2*N*(N-1)
    Elem_idx = Mesh.Edge2Elem(Mesh.DDMData(1,:)==kk,:);
    for j = 1:1/(N/n)
        perm_harm_aver((kk-1)/(N/n)+j) = sum(lambda_S(Elem_idx(j,:)).*coef(Elem_idx(j,:)))/(2*lambda_S(Elem_idx(j,1))*coef(Elem_idx(j,1))*lambda_S(Elem_idx(j,2))*coef(Elem_idx(j,2)));
    end
end