function R = R_covarianceExponential(x1,y1,x2,y2,Lx,Ly,sigma2)

R = sigma2*exp(-abs(x1-x2)/Lx-abs(y1-y2)/Ly);
