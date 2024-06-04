/*
    Computational verifications for Example ??
*/
load "../ZNr-equations.m";

QQ := Rationals();
QQt<t> := FunctionField(QQ);
A3<u,v,z> := AffineSpace(QQ, 3);

// The surface Z(14, 1)
f := ZNrEquations(14,1 : vars:=[u,v,z]);
Z14_1 := Surface(A3, f);

// A curve
A1<T> := AffineSpace(QQ, 1);
C := Curve(Z14_1, [-u*v + u + v^2 - v + 3]);
C := Curve(IrreducibleComponents(C)[1]); //the two components are the plus/minus z
assert Genus(C) eq 0;

// parametrise, then check the families match
aa := eval Read("../Z14-r/infinitefamilies/14-1.m");
E1 := EllipticCurve(aa[1]);
E2 := EllipticCurve(aa[2]);

phi := [-(3*t^2 + 2*t + 4)/(t*(t + 2)),
    -2/t,
    -16*(t + 1)*(3*t^3 - 5*t^2 - 5*t - 11)/(t^2*(t + 2)^2)
];

assert IsInvertible(map<A1 -> C | [Evaluate(p, T) : p in phi]>);

jj, jj_1728 := WNrModuli(14,1 : vars:=phi[1..2]);

assert jj eq jInvariant(E1)*jInvariant(E2);
assert jj_1728 eq (jInvariant(E1) - 1728)*(jInvariant(E2) - 1728);


// TODO: CHECK QUAD TWIST
