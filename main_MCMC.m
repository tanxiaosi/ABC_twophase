
clear


% KL_field
load('KLEeigs.mat')
rng(5);

N = 8;                                 % Number of coarse grids in each direction 
n = 64;                                 % Number of fine grids in each direction
h = 1/n;
mu_w = 1;                               % Viscosity of the water phase
mu_o = 5;                               % Viscosity of the oil phase
add_off = 7;                           % Number of offline basis in each coarse neighborhood
add_on = 3;                            % Number of online basis in each coarse neighborhood 

% perm1 = @(x,varargin)kappa2(x);          % Permeability field 1
% perm2 = @(x,varargin)kappa4(x);          % Permeability field 2

sn = 3;
Nos = sn^trunc;                      % # of snapshots in 2D space

mu = constructsnapshot(sn,trunc);             %create coef for snapshots


S2 = zeros(n^2,1);                    % Initial water saturation
dt = 1;                                 % Time step size
Nt = ceil(1000/dt);                       % Number of time steps

disp(['Offline basis per coarse neighborhood = ',num2str(add_off)])
disp(['Online basis per coarse neighborhood = ',num2str(add_on)])

% W=rand(trunc,1);
% F=Teig*W; 
% coeff = exp(F);

%====================
% initialization
%====================

QuadRule_1D = gauleg(0,1,3);
QuadRule_2D = TProd(gauleg(0,1,6));
O_Handle = @(x,varargin)0*ones(size(x,1),1);
l_Handle = @(x,varargin)ones(size(x,1),1);

clear Mesh

% Initialize function handles for source/sink and boundary data

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

% Assembling fine scale matrices and load vector

[B,mark] = assemMat_Inn2(Mesh,h);
Dofs = Edge2Dof(Mesh,mark);

g = assemLoad_bnd2(Mesh,h,UD1,UD2);
L = assemLoad_Vol_scalar(Mesh,h,QuadRule_2D,G);

M_u = cell(1,trunc);
for i = 1:trunc
    M_u{:,i} = assemMat_Vol(Mesh,h,mark,Teig(:,i),S2,mu_w,mu_o);
end

% Compute basis for the snapshot space

R_snap_u_i = cell(Nos,1);
coef_tot = 0;
M_u_tot = 0;
R_snap_u = [];

for i = 1:Nos
    M_u0 = 0; 
    coef0 = 0;
    for j = 1:trunc
        coef0 = coef0 + Teig(:,j)*mu(i,j);
        M_u0 = M_u0 + mu(i,j)*M_u{:,j}; %size!!!
    end
    coef_tot = coef_tot + coef0;
    M_u_tot = M_u_tot + M_u0;
    [R_snap_u_i{i},R_snap_p] = snapshot_basis(Mesh,M_u0,B,h,Dofs,N,gp);
    R_snap_u = [R_snap_u; R_snap_u_i{i}];
end

for i = 1:Nos
     R_snap_u = [R_snap_u; R_snap_u_i{i}];
end

% POD and Offline space

coef = coef_tot/Nos;

M_ut = M_u_tot/Nos;

temp=[];
for ii = 1:Nos
    for jj = 1:Nos
        temp(:,:,ii,jj) = R_snap_u_i{ii}*M_ut*R_snap_u_i{jj}';
    end
end

for ii = 1:Nos
    for jj = 1:Nos
        newtemp(:,:,ii,jj) = (temp(:,:,ii,jj) + temp(:,:,jj,ii)')/2;
    end
end

M0 = R_snap_u*M_ut*R_snap_u';
M0 = (M0+M0')/2;
perm_harm_aver = coef_aver(Mesh,N,h,coef,S2,mu_w,mu_o);
R_off_u = POD_offline(M0,h,N,add_off,perm_harm_aver);
R_off_u1 =  POD_offline_dis(newtemp,h,N,add_off,perm_harm_aver);

% POD and Online space

rng(5);
theta = randn(trunc, 1);
LogK = Teig*theta;
coef = exp(LogK);

M_u = assemMat_Vol(Mesh,h,mark,coef,S2,mu_w,mu_o);


M0 = R_off_u*R_snap_u*M_u*R_snap_u'*R_off_u';
M0 = (M0+M0')/2;
perm_harm_aver = coef_aver(Mesh,N,h,coef,S2,mu_w,mu_o);
R_on_u = POD_online(M0,h,N,add_on,add_off,perm_harm_aver,R_off_u);

B0 = R_snap_p*B*R_snap_u'*R_off_u'*R_on_u';
A0 = [R_on_u*M0*R_on_u' -B0'; -B0 zeros(size(R_snap_p,1))];
L0 = [zeros(size(R_on_u,1),1);R_snap_p*(g-L)];
aver = zeros(1,size(A0,1));
aver(size(R_on_u,1)+1:end) = 1;
U0 = [A0;aver]\[L0;0];

U_aux = R_snap_u'*R_off_u'*R_on_u'*U0(1:size(R_on_u,1));
p_aux = R_snap_p'*U0(size(R_on_u,1)+1:end);

U_aux = postprocess(Mesh,M_u,B,Dofs,N,gp,U_aux,g-L);

U_on = nodal_value(Mesh,U_aux,mark,h,UD1,UD2);
p_on = nodal_value_scalar(Mesh,p_aux,h);

sat_eq_solver_on
