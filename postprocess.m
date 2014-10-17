function U_new = postprocess(Mesh,M_u,B,Dofs,N,gp,U_aux,L)

U_new = U_aux;

[pnum,primal]=Effect_loc_pu(N,N);

i = N;

Dofsi = Dofs(Mesh.DDMData(1,:)==2*N*(N-1)+i);
Dofsi = Dofsi(Dofsi~=0);
Dofsi_p = find(sum(abs(B(:,Dofsi)'))~=0)';
A_aux = [M_u(Dofsi,Dofsi) -B(Dofsi_p,Dofsi)';-B(Dofsi_p,Dofsi) zeros(numel(Dofsi_p))];
aver = zeros(1,size(A_aux,1));
aver(numel(Dofsi)+1:end) = 1;
Dofs_bd = [];
for k = 1:pnum(i,1)
    pk = primal(i,k);
    Dofs_bd = [Dofs_bd;Dofs(Mesh.DDMData(1,:)==gp(i,pk))];
end
g = [zeros(numel(Dofsi),1); L(Dofsi_p)]-[M_u(Dofsi,Dofs_bd) -B(Dofsi_p,Dofsi)';-B(Dofsi_p,Dofs_bd) zeros(numel(Dofsi_p))]*[U_aux(Dofs_bd);zeros(numel(Dofsi_p),1)];
temp = [A_aux;aver]\[g;0];

U_new(Dofsi) = temp(1:numel(Dofsi));


i = N*(N-1)+1;

Dofsi = Dofs(Mesh.DDMData(1,:)==2*N*(N-1)+i);
Dofsi = Dofsi(Dofsi~=0);
Dofsi_p = find(sum(abs(B(:,Dofsi)'))~=0)';
A_aux = [M_u(Dofsi,Dofsi) -B(Dofsi_p,Dofsi)';-B(Dofsi_p,Dofsi) zeros(numel(Dofsi_p))];
aver = zeros(1,size(A_aux,1));
aver(numel(Dofsi)+1:end) = 1;
Dofs_bd = [];
for k = 1:pnum(i,1)
    pk = primal(i,k);
    Dofs_bd = [Dofs_bd;Dofs(Mesh.DDMData(1,:)==gp(i,pk))];
end
g = [zeros(numel(Dofsi),1); L(Dofsi_p)]-[M_u(Dofsi,Dofs_bd) -B(Dofsi_p,Dofsi)';-B(Dofsi_p,Dofs_bd) zeros(numel(Dofsi_p))]*[U_aux(Dofs_bd);zeros(numel(Dofsi_p),1)];
temp = [A_aux;aver]\[g;0];

U_new(Dofsi) = temp(1:numel(Dofsi));