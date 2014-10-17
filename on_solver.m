% Offline solution

M_u = assemMat_Vol(Mesh,h,mark,coef,S2,mu_w,mu_o);
M0 = R_off_u*R_snap_u*M_u*R_snap_u'*R_off_u';
M0 = (M0+M0')/2;
perm_harm_aver = coef_aver(Mesh,N,h,coef,S2,mu_w,mu_o);
R_on_u = POD_online(M0,h,N,add_on,add_off,perm_harm_aver,R_off_u);

B0 = R_snap_p*B*R_snap_u'*R_off_u'*R_on_u';
A0 = [R_on_u*M0*R_on_u' -B0'; -B0 zeros(size(R_snap_p,1))];
aver = zeros(1,size(A0,1));
aver(size(R_on_u,1)+1:end) = 1;
U0 = [A0;aver]\[L0;0];


U_aux = R_snap_u'*R_off_u'*R_on_u'*U0(1:size(R_on_u,1));
% p_aux = R_snap_p'*U0(size(R_off_u,1)+1:end);

U_aux = postprocess(Mesh,M_u,B,Dofs,N,gp,U_aux,g-L);

% Solve the linear system

% U_off = nodal_value(Mesh,U_aux,mark,h,UD1,UD2);
% p_off = nodal_value_scalar(Mesh,p_aux,h);

% Compute errors

% err1 = L2Err_DGBFE_K(Mesh,U_fine(1:2:end),QuadRule_2D,UD1,coef);
% err2 = L2Err_DGBFE_K(Mesh,U_fine(2:2:end),QuadRule_2D,UD2,coef);
% erru = sqrt(err1^2+err2^2);
% err1 = L2Err_DGBFE_K(Mesh,U_off(1:2:end)-U_fine(1:2:end),QuadRule_2D,UD1,coef);
% err2 = L2Err_DGBFE_K(Mesh,U_off(2:2:end)-U_fine(2:2:end),QuadRule_2D,UD2,coef);
% erru = sqrt(err1^2+err2^2)/erru;
% errp = L2Err_DGBFE(Mesh,p_fine,QuadRule_2D,Pressure);
% errp = L2Err_DGBFE(Mesh,p_off-p_fine,QuadRule_2D,Pressure)/errp;
% disp(['Relative L2 error for velocity (compare with fine solution) is ',num2str(erru)])
% disp(['Relative L2 error for pressure (compare with fine solution) is ',num2str(errp)])