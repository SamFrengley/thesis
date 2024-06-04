/*
    Computational verifications for Example 6.3.
*/
load "../Z12-r/Z12equations.m";
load "../Z2_3_4/iscongruent.m";
QQ := Rationals();
QQt<t> := FunctionField(QQ);
A3<u,v,z> := AffineSpace(QQ, 3);
_<x> := PolynomialRing(QQt);

//The surface Z(12, 1)
Z12_1 := Surface(A3, z^2 - Z12Equations(1, u,v));

//The fibration
Z12_1_t := HyperellipticCurve(Z12Equations(1, x, t));
O := Z12_1_t![0, t^3 - 3*t, 1];
P := Z12_1_t![0, -(t^3 - 3*t), 1];

Z_test := EllipticCurve(Z12_1_t, O);
Z_elliptic := EllipticCurve([-27*(t^8 + 46*t^6 + 859*t^4 - 186*t^2 + 9),
    -54*(t^12 + 69*t^10 + 2082*t^8 + 24731*t^6 - 7848*t^4 + 621*t^2 + 27)]);

assert IsIsomorphic(Z_test, Z_elliptic);

P2t<X,Y,Z> := AmbientSpace(Z_elliptic);
P := Points(Scheme(Z_elliptic, [X - 3*(2*t^8 + 37*t^6 + 450*t^4 - 27*t^2 + 54)/(t^2-3)^2*Z]))[1];
P := Z_elliptic!P;
assert Order(P) eq 0;

//A nicer curve
A1<T> := AffineSpace(QQ, 1);
C := Curve(Z12_1, [2*u + v^2 + 3*v]);
assert Genus(C) eq 0;

aa := eval Read("../Z12-r/infinitefamilies/12-1.m");
E1 := EllipticCurve(aa[1]);
E2 := EllipticCurve(aa[2]);

phi := [-3*t*(t - 3/2)*(t + 2/3)/((t^2 + 1)^(2)),
        -3*t*(t + 2/3)/((t^2 + 1)),
        18*t*(t + 2/3)*(t^2 - 3*t - 1)*(t^2 + 1/2*t + 1/4)/((t^2 + 1)^(3))
        ];

assert IsInvertible(map<A1 -> C | [Evaluate(p, T) : p in phi]>);

jj, jj_1728 := W12Moduli(1, phi[1], phi[2]);
assert jj eq jInvariant(E1)*jInvariant(E2);
assert jj_1728 eq (jInvariant(E1) - 1728)*(jInvariant(E2) - 1728);

assert IsNrCongruent(3, 1, E1, E2);
assert IsNrCongruent(4, 1, E1, E2);
