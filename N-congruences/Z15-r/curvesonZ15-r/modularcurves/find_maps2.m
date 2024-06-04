Attach("../../../ZNr-equations.m");
/* prime := NextPrime(10^30); */
/* F := GF(prime); */
F := Rationals();
A2<u,v> := AffineSpace(F, 2);
P<[x]> := PolynomialRing(Rationals(), 8);


crvs, _ := eval Read("15-2.m");

k := 1;

i := 5;

crv := crvs[i];
m := crv[1]; m;

X := SmallModularCurve(m);  
KX<[x]> := FunctionField(X);
X_K := BaseChange(X, KX);

A2<u,v> := AffineSpace(KX, 2);
jj, jj_1728 := WNrModuli(15,2 : vars:=[u,v]);

C := Curve(A2, Evaluate(crv[k+1], [u,v]));
j_map := map<C -> A2 | [jj,jj_1728]>;

j1 := jInvariant(X_K!(x cat [1]), m);
j2 := jNInvariant(X_K!(x cat [1]), m);


P := A2![j1*j2, (j1-1728)*(j2-1728)] @@ j_map;
