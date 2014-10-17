function [Eigs, W] = getKLEgivenCorrLeng(Lx,Ly,sigma2,Px,Py,numEigs,Ngrid)


% x-direction discretization
dx = Lx/Px;
% y-direction discretization
dy = Ly/Py;

%number of discrete points in the x-direction
Nx = floor(1.0/dx)+1;
%number of discrete points in the x-direction
Ny = floor(1.0/dy)+1;

%fprintf('Nx = %d, Ny = %d \n',Nx,Ny);
N = Nx*Ny;
if(Nx > Ngrid)
    fprintf('Warning in KLE code: Nx > Ngrid, which does not result in good interpolation\n');
    fprintf('  We have Nx=%d, while Ngrid = %d\n', Nx, Ngrid);
end
if(Ny > Ngrid)
    fprintf('Warning in KLE code: Ny > Ngrid, which does not result in good interpolation\n');
     fprintf('  We have Ny=%d, while Ngrid = %d\n', Ny, Ngrid);
end
A = zeros(N,N);
x = 0.5*dx:dx:(Nx-.5)*dx;
y = 0.5*dy:dy:(Ny-.5)*dy;


for k=1:Ny
    for l=1:Nx
        m=(k-1)*(Nx)+l;
        
        for i=1:Ny
            for j=1:Nx
                 R = R_covarianceNormal(x(l),y(k),x(j),y(i),Lx,Ly,sigma2)*dx*dy;
 %               R = R_covarianceExponential(x(l),y(k),x(j),y(i),Lx,Ly,sigma2)*dx*dy;
                if(abs(R) <=1e-10)
                    R = 0;
                end
                n = (i-1)*Nx + j;
                A(m,n) = R;
            end
        end
    end
end


N=size(A,1);

[Vfull, Dfull] = eig(A);

[dtemp I] = sort(diag(Dfull),'descend');
if(numEigs > length(dtemp))
fprintf('ERROR: You did not use a grid size which can produce enough eigenvectors!\n');
fprintf('  You have %d eigenvectors but a %d grid size\n',numEigs,length(dtemp));
end
D = diag( dtemp(1:numEigs) );
V = Vfull(:,I(1:numEigs));

Eigs = dtemp(1:numEigs); % output the eigenvalues to see decay.
% figure; plot(Eigs);

for i=1:numEigs
    if(V(i,1)<0)
        V(i,:) = -V(i,:);
    end
end
if(sum(diag(D))/sigma2 < .7)
fprintf('Warning: The eigenvalues sum to only %4.3f. Use more eigenvalues.\n',sum(diag(D)) );
end
V=V*sqrt(N)*sqrt(D);


% This interpolates V to an Ngrid x Ngrid variable, W
dx=1.0/(Ny-1);
dy=1.0/(Nx-1);
x1=0:dx:1.0;
y1=0:dy:1.0;

[X1,Y1]=meshgrid(x1,y1);
dx=1.0/(Ngrid-1);
x2=0:dx:1.0;
y2=x2;
[X2,Y2]=meshgrid(x2,y2);

S = zeros(Ngrid,Ngrid,numEigs);
for i=1:numEigs
    SS=interp2(X1,Y1,reshape(V(:,i),Nx,Ny),X2,Y2);
    S(:,:,i)=SS(:,:);
end
W=reshape(S,Ngrid*Ngrid,numEigs);
for i=1:numEigs
    if(W(1,i)<=0)
        W(:,i) = -W(:,i);
    end
end
clear S;
clear SS;


