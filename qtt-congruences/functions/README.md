This directory contains `Magma` functions to support the calculations.

- `qttcong.m` contains functions to:

   - Determine if subgroups H+, H of GL<sub>2</sub>(ℤ/Nℤ) satisfies the conditions on Lemma 4.4.1. 
   - Given subgroups H1+,H1 of GL<sub>2</sub>(ℤ/N<sub>1</sub>ℤ) and H2+,H2 of GL<sub>2</sub>(ℤ/N<sub>2</sub>ℤ) satisfying the conditions of Lemma 4.1.1 computes the subgroup Hdelta of GL<sub>2</sub>(ℤ/(N<sub>1</sub>*N<sub>2</sub>)ℤ) such that 
   X^+(H1, H2) = X(Hdelta).

- `mwsieve.m` contains code to perform a Mordell--Weil sieve on a double cover C -> E of an elliptic curve E/ℚ. The data required is a Weierstrass equation for E and an element d of ℚ(E) such that extension ℚ(C)/ℚ(E) is given by adjoining √d.

- `twocover.m` contains code for doing 2-cover descent on a genus 2 curve admitting a Richelot isogeny, as described in Section 4.6. 