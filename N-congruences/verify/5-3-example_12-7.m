/*
    Computational verifications for Example 6.5.
*/
load "../Z12-r/Z12equations.m";
load "../Z2_3_4/iscongruent.m";
QQ := Rationals();
QQt<t> := FunctionField(QQ);
A3<u,v,z> := AffineSpace(QQ, 3);

// The surface Z(12, 7)
Z12_7 := Surface(A3, z^2 - Z12Equations(7, u,v));

// A curve
A1<T> := AffineSpace(QQ, 1);
C := Curve(Z12_7, [u^2 + 3*u - 2*v^2 + 2]);
assert Genus(C) eq 0;

// parametrise, then check the families match
aa := eval Read("../Z12-r/infinitefamilies/12-7.m");
E1 := EllipticCurve(aa[1]);
E2 := EllipticCurve(aa[2]);

phi := [-(t - 3)*(t - 1)*(t + 1)*(t + 3)/(t^4 - 2*t^2 + 9),
        -2*t*(t^2 + 3)/(t^4 - 2*t^2 + 9),
        -2*(t - 3)*(t - 1)*(t + 1)*(t + 3)*(t^2 - 3)*(t^2 + 3)*(t^4 + 3)*(t^4 + 27)/(t^4 - 2*t^2 + 9)^4
        ];

assert IsInvertible(map<A1 -> C | [Evaluate(p, T) : p in phi]>);

jj, jj_1728 := W12Moduli(7, phi[1], phi[2]);
assert jj eq jInvariant(E1)*jInvariant(E2);
assert jj_1728 eq (jInvariant(E1) - 1728)*(jInvariant(E2) - 1728);

assert IsNrCongruent(3, 1, E1, E2);
assert IsNrCongruent(4, 3, E1, E2);