%  This code was written by Tomohiko Mizutani (April 6, 2023)

function [A, W_true, K_true, kappa, omega, beta, condNum] ...
	= genIllcondData(d, n, r, delta, cVal, seed)
    

rng(seed, 'twister');
    
beta = 10^(-cVal / (r-1)); 
S = beta.^(0:r-1);
    
    
W = rand(d, r); 
[F,~,G] = svds(W,r);
    
W = F*diag(S)*G';
W = max(W, 0);
W = getNormalMat(W);
    

kappa = getKappa(W);
omega = getOmega(W);
condNum = getCond(W);



alpha = rand(r,1);
barH = dirRnd(alpha, n-r);
H = [eye(r, r),  barH];

beta = max(barH(:));


R = randn(d, n);
normR = norm(R, 1);

N = (delta/ normR) * R;


V = W * H;
A = V + N;


perm = randperm(n);
%perm =[1:n];
A = A(:, perm);
V = V(:, perm);

[~, idxLst] = sort(perm);
K = idxLst(1:r);

W_true = W;
K_true = K;



end


%%%%%%%%
function [P] = dirRnd(alpha, n)

r = size(alpha, 1);
P = zeros(r, n);

for i=1:r
    
    P(i, :) = gamrnd(alpha(i), 1, 1, n);

end

sumLst = sum(P);
P = P ./ repmat(sumLst, r, 1);

end 


%%%%%%%%
function [A] = getNormalMat(A)
    
[d,n] = size(A);    

isVec = (1 ./ sum(A))';
A = A * spdiags(isVec, 0, n, n);

end

%%%%%%%%
function [condNum] = getCond(W)

    [F, S, G] = svd(W, 'eco');
    condNum = max(diag(S)) / min(diag(S));
    
end


%%%%%%%%
function [kappa] = getKappa(W)
    
[d, r] = size(W);
    
kappa = inf;
for i=1:r
    
    w = W(:,i);
    S = W;
    S(:,i) = [];

    fval = solveLP(S, w);
    
    if fval < kappa
	kappa = fval;
    end
    
end


end


%%%%%%%%
function [fval] = solveLP(A, b)
    
[m,n] = size(A);

% Objective function
f = [zeros(n,1); zeros(m,1); ones(m,1)];

% Inequality constraints: yi - zi <= 0, i = 1, ..., n
C1 = [zeros(m,n), eye(m), -1*eye(m)];
u1 = zeros(m,1);

% Inequality constraints: -yi - zi <= 0, i = 1, ..., n
C2 = [zeros(m,n), -1*eye(m), -1*eye(m)];
u2 = zeros(m,1);

Aineq = [C1; C2];
bineq = [u1; u2];

% Equality constraint: Ax - y = b
Aeq = [A, -1*eye(m), zeros(m,m)];
beq = b;

% Lower and upper bounds
lb = [zeros(n,1); -1*inf*ones(m,1); -1*inf*ones(m,1)];
ub = inf*ones(n+2*m, 1);

% Solve the LP
opt = cplexoptimset('cplex');
opt.display = 'off';
      

[sol, fval] = cplexlp(f, Aineq, bineq, Aeq, beq, lb, ub, [], opt);


    
end

%%%%%%%%
function [omega] = getOmega(W)
    
[d, r] = size(W);

omega = inf;
for i=1:r
   
    col = W(:,i);
    
    S = W;
    S(:,i) = [];
    
    val = min(sum(abs(S - repmat(col, 1, r-1))));

    if val < omega
	omega = val;
    end
    
    
end
    

end






    
    
    











    
    
    









    
    
    








