function gp=Glo2Loc(Nx_subdom,Ny_subdom)

gp = zeros(Nx_subdom*Ny_subdom,4);

for i = 1:Nx_subdom*Ny_subdom
      xy = [ceil(i/Ny_subdom),mod(i,Ny_subdom)];
      if xy(2)==0
          xy(2)=Ny_subdom;
      end
      gp(i,1) = (xy(1)-1)*(Ny_subdom-1) + xy(2)-1;
      gp(i,2) = (xy(1)-1)*(Ny_subdom-1) + xy(2);
      gp(i,3) = Nx_subdom*(Ny_subdom-1) + (xy(1)-2)*Ny_subdom + xy(2);
      gp(i,4) = Nx_subdom*(Ny_subdom-1) + (xy(1)-1)*Ny_subdom + xy(2);
end