function [mu] = constructsnapshot(sn,dim)

N = sn^dim;
mu = zeros(N,dim);

for kk = 1:dim
    t = 3^(dim - kk);
    v = repmat(1:3,[t 1]);
    v = v(:)';
    ind = 3^(kk-1);
    for jj = 1:ind
        a = (jj-1)*N/ind + 1;
        b = jj*N/ind;
        mu(a:b,kk) = v;
    end
end
mu = mu*(1/(sn+1));
    