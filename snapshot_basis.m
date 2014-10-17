function [R_snap_u,R_snap_p] = snapshot_basis(Mesh,M_u,B,h,Dofs,N,gp)

n = 1/h;
Dofsi = Dofs(Mesh.DDMData(1,:)==2*N*(N-1)+1);
Dofsi = Dofsi(Dofsi~=0);
Dofsi_p = find(sum(abs(B(:,Dofsi)'))~=0)';

u_length = numel(Dofsi)+1;
p_length = numel(Dofsi_p);

R_snap_u = zeros(N^2*4*u_length,1);
I_snap_u = R_snap_u;
J_snap_u = R_snap_u;

R_snap_p = ones(N^2*p_length,1);
I_snap_p = R_snap_p;
J_snap_p = R_snap_p;

[pnum,primal]=Effect_loc_pu(N,N);

last_u = 0;
last_p = 0;

for i = 1:N^2
    if mod(i,10)==0
        display(['  processing the ',num2str(i),'-th subdomain'])
    end
    Dofsi = Dofs(Mesh.DDMData(1,:)==2*N*(N-1)+i);
    Dofsi = Dofsi(Dofsi~=0);
    Dofsi_p = find(sum(abs(B(:,Dofsi)'))~=0)';
    A_aux = [M_u(Dofsi,Dofsi) -B(Dofsi_p,Dofsi)';-B(Dofsi_p,Dofsi) zeros(numel(Dofsi_p))];
    aver = zeros(1,size(A_aux,1));
    aver(numel(Dofsi)+1:end) = 1;
    for k = 1:pnum(i,1)
        pk = primal(i,k);
        Dofs_bd = Dofs(Mesh.DDMData(1,:)==gp(i,pk));
        for l = 1:numel(Dofs_bd)
        phi = zeros(size(M_u,1)+size(B,1),1);
        temp = zeros(numel(Dofs_bd)+numel(Dofsi_p),1);
        temp(l) = 1;
        rhs = ones(numel(Dofsi_p),1);
        if pk ==1 || pk ==3
            g = [zeros(numel(Dofsi),1);N^2*h*(2*h)^2*rhs]-[M_u(Dofsi,Dofs_bd) -B(Dofsi_p,Dofsi)';-B(Dofsi_p,Dofs_bd) zeros(numel(Dofsi_p))]*temp;
        else
            g = [zeros(numel(Dofsi),1);-N^2*h*(2*h)^2*rhs]-[M_u(Dofsi,Dofs_bd) -B(Dofsi_p,Dofsi)';-B(Dofsi_p,Dofs_bd) zeros(numel(Dofsi_p))]*temp;
        end
        phi([Dofsi; size(M_u,1)+Dofsi_p]) = [A_aux;aver]\[g;0];
        phi(Dofs_bd(l)) = .5;
        I_snap_u(last_u+1:last_u+u_length) = ((gp(i,pk)-1)*n/N+l)*ones(u_length,1);
        J_snap_u(last_u+1:last_u+u_length) = [Dofsi;Dofs_bd(l)];
        R_snap_u(last_u+1:last_u+u_length) = phi([Dofsi;Dofs_bd(l)]);
        last_u = last_u+u_length;
        end
    end
    I_snap_p(last_p+1:last_p+p_length) = i*ones(p_length,1);
    J_snap_p(last_p+1:last_p+p_length) = Dofsi_p;
    last_p = last_p+p_length;
end

R_snap_u = sparse(I_snap_u,J_snap_u,R_snap_u);
R_snap_p = sparse(I_snap_p,J_snap_p,R_snap_p);

