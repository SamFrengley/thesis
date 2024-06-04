The script `Z14equations.m` records the equations in Theorem ?? together with rational functions jj, (j-1728)(j-1728) \\in Q(u,v) which give the moduli interpretation of W(14, r). A sample usage, which shows that the elliptic curves ?? and ?? are (14,1)-congruent (up to quadratic twist) is as follows:
```
> load "Z14equations.m";
> A3<u,v,z> := AffineSpace(Rationals(), 3);
> f := z^2 - Z14Equations(1, u, v);
> Z14_1 := Surface(A3, f);
...
```

In addition:
- The directory `curvesonZ14-r` contains the data of ?? (in the subdirectory `modularcurves`) and ?? (in the subdirectory `othercurves`).
  
- The directory `infinitefamilies` contains Weierstrass equations for the elliptic curves in ??.
  
- The directory `otherpoints` contains: 
  - In the subdirectory `points`: points we found on the surfaces Z(14, r) (where r \\neq 1) which are not explained by the curves in ??. 
  - In the subdirectory `examples`: the corresponding pairs of (14, r)-congruent elliptic curves (in the same order as the points).