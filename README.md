# RHHP - Refinement of Hottopixx with Hybrid Postprocessing
This is a MATLAB code for the algorithm RHHP, presented in the paper, 

Tomohiko Mizutani, Refinement of Hottopixx Method for Nonnegative Matrix Factorization Under Noisy Separability, SIAM Journal on Matrix Analysis and Applications, 43(3):1029-1057, 2022.

## Requirement
You need to install CPLEX, and include the path of the directiry in which you placed it. See the file 'run_rhhp.m'.

## Quick Start
The following command runs RHHP for a noisy separable matrix  $A = W H + N$ of size d $\times$ n with factorization rank r.


```bash
run_rhhp
```

### Input 
- ``d`` : Number of columns  
- ``n`` : Number of rows
- ``r`` : Factorization rank
- ``delta`` : Noise intensity, i.e, delta = $||N||_1$
- ``dataType`` : Choose 'normal' or 'ill-conditioned'
- ``cVal`` : Parameter for generating $\alpha$ such that 
                     $\alpha^{(\text{r}-1)} = 10^{(-\text{cVal})}$
- ``seed`` : Seed of random numbers
- ``flag_dispMatPara`` : Choose 1 if you display kappa, omega, beta, and condNum (condition number), of input matrix; otherwise, 0

For the details of cVal, kappa, omega and beta, see page 1048 of the paper for cVal; and pages 1031-1032 for kappa, omega and beta.

### Output
- Recovery rate
- Elapsed time

---
Contact: Tomohiko Mizutani [(mizutani.t@shizuoka.ac.jp)](mailto:mizutani.t@shizuoka.ac.jp)
