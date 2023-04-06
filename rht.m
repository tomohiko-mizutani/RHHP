%  Tomohiko Mizutani (April 6, 2023)

function [K, ptLst] = rht(A, r)

    
[d, n] = size(A);    
numVars = n^2 + d*n + 1;

% Objective function
f = sparse(numVars, 1);
f(numVars, 1) = 1;

% Inequality constraint 1: -AX - Y <= -A
C1 = sparse([kron(speye(n,n), -1*A), -1*speye(d*n, d*n), sparse(d*n,1)]);
u1 = sparse(-1*reshape(A, [d*n, 1]));


% Inequalty constraint 2: AX - Y <= A
C2 = sparse([kron(speye(n,n), A), -1*speye(d*n, d*n), sparse(d*n,1)]);
u2 = sparse(reshape(A, [d*n, 1]));


% Inequality constraint 3: Y(1,j) + ... + Y(d,j) - z <= 0, j = 1, ..., n
C3 = sparse([sparse(n, n^2), kron(speye(n,n), ones(1,d)), -1*ones(n,1)]);
u3 = sparse(n,1);


% Inequality constraint 4: -X(i,i) + X(i,j) <= 0, i,j = 1, ..., n
idxLst = [1:n^2];
idxLst(1, [1:n+1:n^2]) = -1;

L = reshape(idxLst, [n,n])';
I = reshape(L, [n^2, 1]);
I(I < 0) = [];

Cpos = [sparse([1:n^2-n], I, 1, n^2-n, n^2), sparse(n^2-n, d*n+1)];

I = kron([1:n+1:n^2], ones(1,n-1));
Cneg = [sparse([1:n^2-n], I, -1, n^2-n, n^2), sparse(n^2-n, d*n+1)];

C4 = Cpos + Cneg;
u4 = sparse(n^2-n, 1);

Aineq = [C1;C2;C3;C4];
bineq = [u1;u2;u3;u4];



% Equality constraint: X(1,1) + ... + X(n,n) = r
I = [1:n+1:n^2];
Aeq = [sparse(1, I, 1, 1, n^2), sparse(1, d*n+1)];
beq = [r];

% Lower bound on variable
lb = sparse([sparse(n^2, 1); -1*inf*ones(d*n+1, 1)]);


% Upper bound on variable
ub = sparse(inf*ones(numVars,1));
I = [1:n+1:n^2];
ub(I,1) = 1;



% Solve the model by cplex
opt = cplexoptimset('cplex');
opt.display = 'off';
%opt.Algorithm = 'interior-point';


sol = cplexlp(f, Aineq, bineq, Aeq, beq, lb, ub, [], opt);


x = sol([1:n^2],1);
xMat = reshape(x, [n,n]);
ptLst = diag(xMat);


[~, idxLst] = sort(ptLst, 'descend');
K = idxLst([1:r])';

    
end