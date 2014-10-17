% Run script for two phase-flow simulation (water and oil)

% Initialize constants
tic
N = 8;                                 % Number of coarse grids in each direction 
n = 64;                                 % Number of fine grids in each direction
h = 1/n;
mu_w = 1;                               % Viscosity of the water phase
mu_o = 5;                               % Viscosity of the oil phase
add_off = 7;                           % Number of offline basis in each coarse neighborhood
add_on = 3;                            % Number of online basis in each coarse neighborhood 

perm1 = @(x,varargin)kappa2(x);          % Permeability field 1
perm2 = @(x,varargin)kappa4(x);          % Permeability field 2
mu = [.3 .4 .6 .9];
mu_real = .5;

S2 = zeros(n^2,1);                    % Initial water saturation
dt = 1;                                 % Time step size
Nt = ceil(1000/dt);                       % Number of time steps
switchf = 0;                               % switch for fine solver 0 off 1 on
switchm = 1;                               % switch for ms solver

disp(['Offline basis per coarse neighborhood = ',num2str(add_off)])
disp(['Online basis per coarse neighborhood = ',num2str(add_on)])

% Run the fine scale solver for velocity and saturation

if switchf == 1
    fine_solver
    sat_eq_solver_fine
end

% return;
% Run the multiscale solver for velocity and saturation

if switchm == 1
    snap_off_on_solver
    sat_eq_solver_on
end

% Compute the relative L2 error of the saturation at final time
return;

if para == 1
    errS = L2Err_DGBFE(Mesh,S_on-S_fine,QuadRule_2D,O_Handle)/...
        L2Err_DGBFE(Mesh,S_fine,QuadRule_2D,O_Handle);
    disp(['Relative L2 error for saturation at final time is '...
        ,num2str(errS)]);
end
toc