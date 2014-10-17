% Run script for fine scale biconstant finite element solver.

% Chak Shing Lee, October 2013
% Department of Mathematics
% Texas A&M University

% clear all

QuadRule_1D = gaulob(0,1,3);
QuadRule_2D = TProd(gauleg(0,1,6));
O_Handle = @(x,varargin)0*ones(size(x,1),1);
l_Handle = @(x,varargin)ones(size(x,1),1);

clear Mesh

G = @(x,varargin)ones(size(x,1),1).*(x(:,2)>(1-2*h)).*(x(:,1)<2*h)-ones(size(x,1),1).*(x(:,2)<2*h).*(x(:,1)>(1-2*h));
UD1 = O_Handle;
UD2 = O_Handle;
Pressure = O_Handle;

% Initialize mesh

Mesh = TProd_Mesh(0:h:1);

Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh = add_Edge2Elem(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;
Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
gp=Glo2Loc(N,N);
Mesh = add_DDMData(Mesh,1,1,N,N,gp);

perm = @(x,varargin)mu_real*perm1(x)+(1-mu_real)*perm2(x);
%coef = perm_k(perm,Mesh);
%coef = reshape(exp(LogK),4096,1);
%coef = perm_k(perm,Mesh);
%coef = exp(LogK);
% [coef,~] = copy_sb10(Mesh,N);

% Assembling matrices and load vector

[B,mark] = assemMat_Inn2(Mesh,h);
Dofs = Edge2Dof(Mesh,mark);
M_u = assemMat_Vol(Mesh,h,mark,coef,S2,mu_w,mu_o);
g = assemLoad_bnd2(Mesh,h,UD1,UD2);

L = assemLoad_Vol_scalar(Mesh,h,QuadRule_2D,G);
L = -g+L;
A = B*assemMat_inv(M_u)*(B');
A = [A;ones(1,size(A,1))];

% Solve the linear system

p_aux = A\[L;0];
U_aux = assemMat_inv(M_u)*(B')*p_aux;


% Plot out solution

%         U_fine = nodal_value(Mesh,U_aux,mark,h,UD1,UD2);
%         p_fine = nodal_value_scalar(Mesh,p_aux,h);
%         plot_BFE(p,Mesh);
%         colorbar;
%         plot_BFE(U(1:2:end),Mesh);
%         colorbar;
%         plot_BFE(U(2:2:end),Mesh);
%         colorbar;

