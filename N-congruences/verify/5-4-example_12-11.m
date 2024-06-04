/*
    Computational verifications for Example 6.6.
*/
load "../Z12-r/Z12equations.m";
load "../Z2_3_4/iscongruent.m";
QQ := Rationals();
QQt<t> := FunctionField(QQ);
A3<u,v,z> := AffineSpace(QQ, 3);

// The surface Z(12, 11)
Z12_11 := Surface(A3, z^2 - Z12Equations(11, u,v));

// A curve
A1<T> := AffineSpace(QQ, 1);
C := Curve(Z12_11, [(9*v^4 + 30*v^2 + 5)*u - 45*v^4 - 6*v^2 + 7]);
C := Curve(IrreducibleComponents(C)[1]); //the two components are the plus/minus z
assert Genus(C) eq 0;

// load the recorded family
aa := eval Read("../Z12-r/infinitefamilies/12-11.m");
E1 := EllipticCurve(aa[1]);
E2 := EllipticCurve(aa[2]);

// parametrise, then check the families match
phi := [-(t^2 - 3)*(7*t^2 + 15)/(5*t^4 + 30*t^2 + 9),
        -1/t,
        -2*(t - 1)*(t + 1)*(t^2 + 3)*(t^4 - 18*t^2 - 27)*(t^4 + 6*t^2 - 3)*(25*t^6 + 63*t^4 + 
        27*t^2 - 27)/((t)^(6)*(5*t^4 + 30*t^2 + 9)^(3))
        ];

assert IsInvertible(map<A1 -> C | [Evaluate(p, T) : p in phi]>);

jj, jj_1728 := W12Moduli(11, phi[1], phi[2]);
assert jj eq jInvariant(E1)*jInvariant(E2);
assert jj_1728 eq (jInvariant(E1) - 1728)*(jInvariant(E2) - 1728);

assert IsNrCongruent(3, 2, E1, E2);
assert IsNrCongruent(4, 3, E1, E2);