This directory contains data to do with 2, 3, and 4-congruences. In particular:

- The script `iscongruent.m` contains a two functions `IsNrCongruent(N, r, E1, E2)` and `IsNrCongruentToTwist(N, r, E1, E2)`. The first employs Theorem 3.3 to check whether a pair of elliptic curves are (N,r)-congruent and the second checks if E<sub>1</sub> and E<sub>2</sub> are congruent up to twist, and if so returns the squarefree d giving the twist(s).
  
- The file `polys.m` records the poynomials from Theorem 3.3 in functions `CongPoly(N,r)` and `TauPoly(N,r)`.
  
- The subdirectory `families` contains files `N-r-congruent.m` which have symmetrical 2-parameter families of elliptic curves E<sub>1</sub>/Q(a,b) and E<sub>2</sub>/Q(a,b) which are (N, r)-congruent for N \\leq 4. The involution swapping E<sub>1</sub> and E<sub>2</sub> is given by swapping a and b. These families are calculated by parametrising the modular diagonal quotient surface Z(N, r).