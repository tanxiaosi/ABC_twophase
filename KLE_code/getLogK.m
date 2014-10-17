clear all
close all


Lx = 1;
Ly = 1;
sigma2 = 2;
Px=5;
Py=5;
mu = 5; % KL-expansion mean.

NumEigs = 5;
Ngrid = 60;


[W] = getKLEgivenCorrLeng(Lx, Ly, sigma2, Px, Py, NumEigs, Ngrid);

theta = randn(NumEigs, 1);

LogK = reshape(mu+ W*theta, Ngrid, Ngrid)';

figure;imagesc(LogK);