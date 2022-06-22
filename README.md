# approximation-of-matrix-norms
Matlab code repository for numerical experiments from the paper "Tight computationally efficient approximation of matrix norms with applications" 
by A. Juditsky, G. Kotsalis and A. Nemirovski:

1. "Boeing" peak-to-peak optimized control experiment from Section 3.3.3 of the paper
    
    Code in /boeing folder; use boeing_demo.m to start the simulation

2. System identification experiment from Section 4.3.2
    
    Code in /ident folder; use id_demo.m to start the simulation
 
To run the experiments Mosek optimization solver [1] and CVX [2] should be installed with appropriate license.

[1] E. D. Anderson. The MOSEK optimization toolbox for MATLAB manual. https://docs.mosek.com/latest/toolbox/index.html

[2] M. Grant and S. Boyd. The CVX Users Guide. http://web.cvxr.com/cvx/doc/
