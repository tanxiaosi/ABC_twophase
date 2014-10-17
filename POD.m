function R_off_u = POD(M0,B0,h,N,add)

R_off_u = zeros(N*(N-1)*2*add,N*(N-1)*2/(N*h));

for kk = 1:2*N*(N-1)
    M_local = M0((kk-1)/(N*h)+1:kk/(N*h),(kk-1)/(N*h)+1:kk/(N*h));
    B_local = B0(:,(kk-1)/(N*h)+1:kk/(N*h));
    X_local = M_local + B_local'*B_local;
    if add < size(M_local,1)
%         [V,lam_i] = eigs(X_local,M_local,add,'LM');
        [V,lam_i]=eig(full(X_local),full(M_local));
        [~,I]=sort(diag(lam_i),'descend');
        V=V(:,I(1:add));
    else
        [V,~] = eig(full(X_local),full(M_local));
    end
%         disp(lam_i)
    R_off_u((kk-1)*add+1:kk*add,(kk-1)/(N*h)+1:kk/(N*h)) = V';
end
