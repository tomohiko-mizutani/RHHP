%  Tomohiko Mizutani (April 6, 2023)

function [K, timeRHHP] = rhhp(A, r)

elpStart = tic;

[K_rht, ptLst] = rht(A, r);
K_rpp = rpp(A, r, ptLst);

resErr1 = compResErr(A, A(:, K_rht), r);

if min(K_rpp) == 0
    resErr2 = inf;
else
    resErr2 = compResErr(A, A(:, K_rpp), r);
end


if resErr1 < resErr2
    K = K_rht;
else
    K = K_rpp;
end
    

timeRHHP = toc(elpStart);


end


%%%%%%%%
function [resErr] = compResErr(A, B, r)

epsilon = 1.0e-8;
    
[d,n] = size(A);    
X = [];

for i=1:n
    
    a = A(:,i);
    sol = cplexlsqnonneglin(B, a);    
    X(:,i) = sol;
    
end

resErr = norm(A-B*X, 'fro')^2;
    
    
end


