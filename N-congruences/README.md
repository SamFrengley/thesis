This directory is effectively a local copy of that attached to my [12-congruences](https://github.com/SamFrengley/12-congruences.git) paper, with 14 and 15 also covered as cases. All of the code is written in [`Magma`](http://magma.maths.usyd.edu.au/magma/).

The most useful resource is `ZNr-equations.m` which provides intrinsics to get the moduli interpretation of the surfaces $Z_{N,r}$. In particular we give equations for $Z_{N,r}$ for each $(N,r)$ in the thesis, as well as recalling those of Fisher $(13,r)$ and $(17,r)$. We also give models in the cases $(18,r)$, $(20,1)$, $(20,3)$, $(22,1)$, $(24,11)$, and $(24,23)$. Note however that **we do not give proofs** except in the cases specified in the thesis. These equations can be accessed with `ZNrEquations` and `WNrModuli`.
  
The directories are arranged as follows:
- In the directory `verify` we provide code which checks the claims made in the paper.

- The directories `Z2_3_4` and `Z5-r` contains files related to 2, 3, 4, and 5-congruences of elliptic curves.

- The directories `Z12-r`, `Z14-r`, `Z15-r` contains files related to $N$-congruences. We record examples of rational points and curves of low genus on these surfaces, and the corresponding pairs of $(N, r)$-congruent elliptic curves.
