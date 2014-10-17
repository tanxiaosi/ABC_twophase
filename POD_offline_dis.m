function R_off_u = POD_offline_dis(newtemp,h,N,add,coef_aver)

n = 1/h;
% R_off_u = zeros(N*(N-1)*2*add,4*N*(N-1)*2/(N/n));
R_off_u = zeros(N*(N-1)*2*add,9*N*(N-1)*2/(N/n));

for kk = 1:2*N*(N-1)
    Index = [(kk-1)/(N/n)+1:kk/(N/n) ((kk-1)/(N/n)+1:kk/(N/n))+2*N*(N-1)/(N/n)...
        ((kk-1)/(N/n)+1:kk/(N/n))+4*N*(N-1)/(N/n) ((kk-1)/(N/n)+1:kk/(N/n))+6*N*(N-1)/(N/n)...
        ];
%     for c = 1:9
%         Index(((c-1)*8+1):c*8) = ((kk-1)/(N/n)+1:kk/(N/n))+2*(c-1)*N*(N-1)/(N/n);
%     end
    Index1 = (kk-1)/(N/n)+1:kk/(N/n);
    temp_store = cell(4,1);
    for ii = 1:4
        for jj = 1:4
            cc = newtemp(:,:,ii,jj);
            temp_store{ii} = [temp_store{ii} cc(Index1,Index1)];
        end
    end
    M_local = [temp_store{1};temp_store{2};temp_store{3};temp_store{4}];
    M_local = sparse(M_local);
%     M_local = M0(Index,Index);
    %     X_local = M0_1((kk-1)/(N*h)+1:kk/(N*h),(kk-1)/(N*h)+1:kk/(N*h))/h;
%     I_local = diag(coef_aver((kk-1)/(N/n)+1:kk/(N/n)));
%     I_local_snap = diag([coef_aver((kk-1)/(N/n)+1:kk/(N/n));coef_aver((kk-1)/(N/n)+1:kk/(N/n));coef_aver((kk-1)/(N/n)+1:kk/(N/n));coef_aver((kk-1)/(N/n)+1:kk/(N/n))])/n;
%     I_local_snap = [I_local I_local I_local I_local;I_local I_local I_local I_local;I_local I_local I_local I_local;I_local I_local I_local I_local];
    %     [V,lam_i] = eigs(X_local,M_local,add,'LM');
    [V,lam_i]=eig(full(M_local));
    [lam_i,I]=sort(diag(lam_i),'descend');
    V=V(:,I(1:add));
    
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
    R_off_u((kk-1)*add+1:kk*add,Index) = V';
    if isreal(R_off_u)~=1
        kk
        break;
    end
end
[I,J,K] = find(R_off_u);
R_off_u = sparse(I,J,K);