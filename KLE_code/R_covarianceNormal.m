function R = R_covarianceNormal(x1,y1,x2,y2,Lx,Ly,sigma2)
R = sigma2*exp(-(x1-x2).^2/(2*Lx.^2)-(y1-y2).^2/(2*Ly.^2));
