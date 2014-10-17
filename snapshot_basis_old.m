function [R_snap_u,R_snap_p] = snapshot_basis(Mesh,M_u,B,h,Dofs,N,gp)

n = 1/h;
R_snap_u = sparse(2*N*(N-1)*(n/N),size(M_u,1));
R_snap_p = sparse(N^2,size(B,1));

[pnum,primal]=Effect_loc_pu(N,N);

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
        if pk ==1 || pk ==3
            g = [zeros(numel(Dofsi),1);N^2*(2*h)^2*h*ones(numel(Dofsi_p),1)]-[M_u(Dofsi,Dofs_bd) -B(Dofsi_p,Dofsi)';-B(Dofsi_p,Dofs_bd) zeros(numel(Dofsi_p))]*temp;
        else
            g = [zeros(numel(Dofsi),1);-N^2*(2*h)^2*h*ones(numel(Dofsi_p),1)]-[M_u(Dofsi,Dofs_bd) -B(Dofsi_p,Dofsi)';-B(Dofsi_p,Dofs_bd) zeros(numel(Dofsi_p))]*temp;
        end
        phi([Dofsi; size(M_u,1)+Dofsi_p]) = [A_aux;aver]\[g;0];
        phi(Dofs_bd(l)) = .5;
        R_snap_u((gp(i,pk)-1)*n/N+l,:) = R_snap_u((gp(i,pk)-1)*n/N+l,:) + phi(1:size(M_u,1))';
        end
    end
    R_snap_p(i,Dofsi_p) = 1;
end

