function M = Local_Solver(Mesh,M_u,M_q,B,h,Dofs,Nx_subdom,Ny_subdom)

% Pre-allocate memory

M = sparse(size(M_u,1),size(M_u,1));

% Local solvers

for i = 1:Nx_subdom*Ny_subdom
    if mod(i,10)==0
        display(['  processing the ',num2str(i),'-th subdomain'])
    end
    Dofse = Dofs(Mesh.DDMData(i,:)==1);
    Dofse = Dofse(Dofse~=0);
    Pj = sparse([1:size(Dofse,1) size(Dofse,1)],[Dofse;size(M_u,1)],[ones(size(Dofse,1),1);0]);
    A_local = M_u(Dofse,Dofse)+B(:,Dofse)'*M_q*B(:,Dofse);
    Aloc = Pj'*(A_local\Pj);
    M = M+Aloc;
end