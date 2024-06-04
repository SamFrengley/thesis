The script `Z12equations.m` records the equations in Theorem 1.1 together with rational functions jj, (j-1728)(j-1728) \\in Q(u,v) which give the moduli interpretation of W(12, r). A sample usage, which shows that the elliptic curves [`735.d2`](https://www.lmfdb.org/EllipticCurve/Q/735/d/2) and [`24255.bh1`](https://www.lmfdb.org/EllipticCurve/Q/24255/bh/1) are (12,7)-congruent (up to quadratic twist) is as follows:
```
> load "Z12equations.m";
> A3<u,v,z> := AffineSpace(Rationals(), 3);
> f := z^2 - Z12Equations(7, u, v);
> Z12_7 := Surface(A3, f);
> P := Z12_7![ -5/13, 21/13, 3402/2197 ];
> jj, jj_1728 := W12Moduli(7, -5/13, 21/13);
> j1 := jInvariant(EllipticCurve("735c1"));
> j2 := jInvariant(EllipticCurve("24255bq2"));
> assert jj eq j1*j2;
> assert jj_1728 eq (j1 - 1728)*(j2 - 1728);
```

In addition:
- The directory `curvesonZ12-r` contains the data of Tables 2--5 (in the subdirectory `modularcurves`) and Tables 6--8 (in the subdirectory `othercurves`).
  
- The directory `infinitefamilies` contains Weierstrass equations for the elliptic curves in Examples 6.3--6.6.
  
- The directory `otherpoints` contains: 
  - In the subdirectory `points`: points we found on the surfaces Z(12, r) (where r \\neq 1) which are not explained by the curves in Tables 2--8. 
  - In the subdirectory `examples`: the corresponding pairs of (12, r)-congruent elliptic curves (in the same order as the points).
  
- The directory `smallexamples` contains the data from Table 9 and a list of (12, 1)-congruent elliptic curves where both curves are contained in the LMFBD (in fact it also has the more refined data of the corresponding points on Z(12, r) and Weierstrass equations for the corresponding elliptic curves).