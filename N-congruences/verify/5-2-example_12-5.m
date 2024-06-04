/*
    Computational verifications for Example 6.4.
*/
load "../Z12-r/Z12equations.m";
load "../Z2_3_4/iscongruent.m";
QQ := Rationals();
calC := EllipticCurve("240c2"); KK<xi,eta> := FunctionField(calC); 
A3<u,v,z> := AffineSpace(QQ, 3);

//The surface Z(12, 5)
Z12_5 := Surface(A3, z^2 - Z12Equations(5, u,v));

//A curve
C := Curve(Z12_5, [u^3 + 3*u^2 - (v^2 + 10)*u - 3*v^2 + 6]);
O := C![-1, 3, 300];
assert Genus(C) eq 1;

// load the recorded family
aa := eval Read("../Z12-r/infinitefamilies/12-5.m");
E1 := EllipticCurve(aa[1]);
E2 := EllipticCurve(aa[2]);

// isomorphism, then check the families match
_, phi1 := EllipticCurve(ProjectiveClosure(C), O);
_, phi2 := IsIsomorphic(Codomain(phi1), calC);

phi := map<C -> calC | [Evaluate(p, [C.1, C.2, C.3, 1]) : p in DefiningPolynomials(phi1*phi2)]>;
_, phi := IsInvertible(phi);
phi := [KK!p : p in DefiningPolynomials(phi)];

jj, jj_1728 := W12Moduli(5, phi[1], phi[2]);
assert jj eq jInvariant(E1)*jInvariant(E2);
assert jj_1728 eq (jInvariant(E1) - 1728)*(jInvariant(E2) - 1728);

assert IsNrCongruent(3, 2, E1, E2);
assert IsNrCongruent(4, 1, E1, E2);