function M1 = assemMat_inv(M)

M1 = 0*M;
for i = 1:size(M,1)
    M1(i,i) = 1/M(i,i);
end