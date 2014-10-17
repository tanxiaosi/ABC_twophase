% Run script for fine scale biconstant finite element solver.

% Chak Shing Lee, October 2013
% Department of Mathematics
% Texas A&M University

M_u = assemMat_inv(assemMat_Vol(Mesh,h,mark,coef,S2,mu_w,mu_o));
A = B*M_u*(B');
A = [A;ones(1,size(A,1))];

% Solve the linear system

p_aux = A\[Lg;0];
U_aux = M_u*(B')*p_aux;
% U_fine = nodal_value(Mesh,U_aux,mark,h,UD1,UD2);
% p_fine = nodal_value_scalar(Mesh,p_aux,h);

% Plot out solution

%         plot_BFE(p,Mesh);
%         colorbar;
%         plot_BFE(U(1:2:end),Mesh);
%         colorbar;
%         plot_BFE(U(2:2:end),Mesh);
%         colorbar;

