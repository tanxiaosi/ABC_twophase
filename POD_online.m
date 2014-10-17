function R_on_u = POD_online(M0,h,N,add_on,add_off,coef_aver,R_off_u)

n = 1/h;
R_on_u = zeros(N*(N-1)*2*add_on,N*(N-1)*2*add_off);

for kk = 1:2*N*(N-1)
    M_local = M0((kk-1)*add_off+1:kk*add_off,(kk-1)*add_off+1:kk*add_off);
    %     X_local = M0_1((kk-1)/(N*h)+1:kk/(N*h),(kk-1)/(N*h)+1:kk/(N*h))/h;
    R_off_local = R_off_u((kk-1)*add_off+1:kk*add_off,(kk-1)/(N/n)+1:kk/(N/n));
    I_local = R_off_local*diag(coef_aver((kk-1)/(N/n)+1:kk/(N/n)))*R_off_local'/n;
    %         [V,lam_i] = eigs(X_local,M_local,add,'LM');
    [V,lam_i]=eig(full(I_local),full(M_local));
    [lam_i,I]=sort(diag(lam_i),'ascend');
    V=V(:,I(1:add_on));
    
%         if kk == 3
%     %         full(X_local)
%             disp(lam_i')
%     %         plot(lam_i)
%     %         hold on
%     %         plot(lam_i,'r*')
%     
%     %         V
%         end
    %     V = [ones(size(V,1),1) V(:,1:end-1)];
    R_on_u((kk-1)*add_on+1:kk*add_on,(kk-1)*add_off+1:kk*add_off) = V';
end
[I,J,K] = find(R_on_u);
R_on_u = sparse(I,J,K);