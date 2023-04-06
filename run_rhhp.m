%  Tomohiko Mizutani (April 6, 2023)

function [] = run_rhhp()

% 
% ### DESCRIPTION ###
% This code runs the algorithm RHHP, presented in the paper,
% T. Mizutani, Refinement of Hottopixx Method for Nonnegative Matrix Factorization
% Under Noisy Separability, SIMAX, INPUT.
%
%    
% ### INPUT ###
% A = W * H + N : noisy separable matrix of size d * n with factorization rank r
%     
% d                 : Number of rows
% n                 : Number of columns
% r                 : Factorization rank
% delta             : Noise intensity, i.e, delta = ||N||_1
% dataType          : Choose 'normal' or 'ill-conditioned'
% cVal              : Parameter for generating alpha such that alpha^(r-1) = 10^(-cVal)
% seed              : Seed of random numbers
% flag_dispMatPara  : Choose 1 if you display kappa, omega, beta
%    
% For the details of cVal, kappa, omega and beta, see page 1048 of the paper for cVal;
% and pages 1031-1032 for kappa, omega and beta.
%    
%
% ### OUTPUT ###
% recovery rate and elapsed time
%    
%  
% ### REQUIREMENT ###
% You need to install CPLEX, and then include the path. See below.
%    


% Add the path of the directory in which you installed CPLEX
addpath(['']);

    
d = 10;
n = 30;
r = 5;
delta = 0.05;

dataType = 'normal';
%dataType = 'ill-conditioned';
%cVal = 4;
seed = 276;

flag_dispMatPara = 0;



switch dataType
    
  case 'normal'

    if flag_dispMatPara == 1
    
    	fprintf('[Data type] \n');
    	fprintf('Normal \n');

    end

    [A, W_true, K_true, kappa, omega, beta, condNum] ...
        = genNormalData(d, n, r, delta, seed);
    
  case 'ill-conditioned'

    if flag_dispMatPara == 1    

    	fprintf('[Data type] \n');
    	fprintf('Ill-conditioned \n');

    end

    [A, W_true, K_true, kappa, omega, beta, condNum] ...
        = genIllcondData(d, n, r, delta, cVal, seed);
    
end

if flag_dispMatPara == 1    

    dispMatPara(kappa, omega, beta, condNum);
    fprintf('\n');

end


[K, timeRHHP] = rhhp(A, r);
rate = compRecoveryRate(K, K_true);


fprintf('[Results] \n');
fprintf('Recovery Rate: %1.3f \n', rate);
fprintf('Elapsed Time: %1.1f s \n', timeRHHP);



end

%%%%%%%%
function rate = compRecoveryRate(tIdxSet, J);

corrNum = size(intersect(tIdxSet, J), 2);
rate = corrNum / size(tIdxSet, 2);

end



%%%%%%%%
function [] = dispMatPara(kappa, omega, beta, condNum);

fprintf('\n');
fprintf('[Parameter values of input matrix] \n');
fprintf('(kappa, omega, beta, condNum) = (%1.4f, %1.4f, %1.4f, %1.4f) \n', ...
	kappa, omega, beta, condNum);


end



