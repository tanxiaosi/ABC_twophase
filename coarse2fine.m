 c2f=sparse(n^2,(n/2)^2);

%     [ix,iy]=[ceil(i/Nf),mod(i,Nf)];
%     [ix,iy]=[ceil(ix/2),ceil(iy/2)];
  
for i = 1:n/2
    for j = 1:n/2
        ref_no = 1+2*(j-1)+2*n*(i-1);
        ind=[ref_no,ref_no+1,ref_no+n,ref_no+n+1];
        c2f(ind,i+(j-1)*n/2)=1;
    end
end
%S= c2f*S;

% S=linspace(0,1,256);
% S2 = nodal_value_scalar(Mesh,S,h);
% S= c2f*S';
% (S2-S)'