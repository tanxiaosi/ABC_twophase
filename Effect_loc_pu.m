function [pnum,primal]=Effect_loc_pu(Nx_subdom,Ny_subdom)


  % Initialize constants
  
  pnum = zeros(Nx_subdom*Ny_subdom,1);
  primal = zeros(Nx_subdom*Ny_subdom,4);
  
  for i = 1:Nx_subdom*Ny_subdom
      xy = [ceil(i/Ny_subdom),mod(i,Ny_subdom)];
      if xy(2)==0
          xy(2)=Ny_subdom;
      end
      if xy(2)~=1
          pnum(i,1) = pnum(i,1)+1;
          primal(i,pnum(i,1)) = 1;
      end
      if xy(2)~=Ny_subdom
          pnum(i,1) = pnum(i,1)+1;
          primal(i,pnum(i,1)) = 2;
      end
      if xy(1)~=1
          pnum(i,1) = pnum(i,1)+1;
          primal(i,pnum(i,1)) = 3;
      end
      if xy(1)~=Nx_subdom
          pnum(i,1) = pnum(i,1)+1;
          primal(i,pnum(i,1)) = 4;
      end
  end