% A two-dimensional Log-Gaussian process. 
%--------------------------------------------------------------------------
close all; clear all;
% tic

addpath ./KLE_code;

Lx = 0.1; % correlation length in x direction.
Ly = 0.1; % correlation length in y direction.
sigma2 = 2;
Px= 2; % Px/Lx: the # of grid mesh before interpolation. 
Py= 1;
mu = 0; % KL-expansion mean.

trunc = 2;
Ngrid = 64;


[Eigs Teig] = getKLEgivenCorrLeng(Lx, Ly, sigma2, Px, Py, trunc, Ngrid);
% figure;plot(Eigs,'-*');

theta = randn(trunc, 1);
LogK = reshape(mu*ones(1,1)+ Teig*theta, Ngrid, Ngrid)';
figure;imagesc(LogK);colorbar;


save('KLEeigs.mat','mu','Teig','trunc','Ngrid');

% time=toc