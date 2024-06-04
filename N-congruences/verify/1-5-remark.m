/*
    Checking the claims in Remark 3.4. 
*/

// Note that curves here have Cremona labels, in the paper they have LMFDB labels
load "../Z2_3_4/polys.m"; 
load "../Z2_3_4/iscongruent.m";

function ProjectiveRoots(f)
    // Given a homogeneous polynomial f \in Q[x,y] return
    // the points in P1 which vanish on f
    _<x> := PolynomialRing(Rationals());
    f1 := Evaluate(f, [x,1]);
    f2 := Evaluate(f,[1,x]);
    ret := [[a[1],1] : a in Roots(f1)] cat [[1,a[1]] : a in Roots(f2)];
    return ret;
end function;

// ----------------
// (i): The case when (N,r) = (3,1)
KK<J1,J2> := FunctionField(Rationals(), 2);
PP<x> := PolynomialRing(KK);
LL<alpha> := quo<PP | x^3 - 3*J1*J2*x - J1*J2*(J1 + J2)>;
alpha_prime := 1/(J1 - J2)*(alpha^2/J2 - alpha - 2*J1);
assert alpha_prime^3 eq J1/J2;

// ----------------
// (ii): The case when (N,r) = (4,1)
PP<x> := PolynomialRing(Rationals());
E1 := EllipticCurve("196b2"); J1 := jInvariant(E1)/1728;
E2 := EllipticCurve("196b1"); J2 := jInvariant(E2)/1728;

Pa,Pb := CongPolys(2,1);
assert exists(alpha){a[1] : a in Roots(Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), x, 0])) | a[1] lt 0};
assert #Roots(Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), alpha, x])) eq 0;

// ----------------
// (ii): The case (N,r) = (3,1)
PP<x> := PolynomialRing(CyclotomicField(3));
E1 := EllipticCurve("15a5"); J1 := jInvariant(E1)/1728;
E2 := EllipticCurve("15a1"); J2 := jInvariant(E2)/1728;

Pa,Pb := CongPolys(3,1);
rts_a := Roots(Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), x, 0]));
assert exists(alpha){a[1] : a in rts_a | #[x : x in Eltseq(a[1]) | x ne 0] gt 1};
assert #Roots(Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), alpha, x])) eq 0;

// ----------------
// (ii): The case (N,r) = (3,2)
PP<x> := PolynomialRing(CyclotomicField(3));
E1 := EllipticCurve("11a2"); J1 := jInvariant(E1)/1728;
E2 := EllipticCurve("11a1"); J2 := jInvariant(E2)/1728;

Pa,Pb := CongPolys(3,2);
rts_a := Roots(Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), x, 0]));
assert exists(alpha){a[1] : a in rts_a | #[x : x in Eltseq(a[1]) | x ne 0] gt 1};
assert #Roots(Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), alpha, x])) eq 0;

// ----------------
// (iii): The case when (N,r) = (3,1)
// (cbrt(J) + 1)(cbrt(J') + 1) case
PP<x> := PolynomialRing(Rationals());
E1 := EllipticCurve("245a1"); J1 := jInvariant(E1)/1728;
E2 := EllipticCurve("1323i1"); J2 := jInvariant(E2)/1728;

Pa,Pb := CongPolys(3,1);
alpha := Roots(Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), x, 0]))[1][1];
assert #Roots(Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), alpha, x])) gt 0; // have alpha, beta
assert GF(3)!(TraceOfFrobenius(E1, 11) - TraceOfFrobenius(E2, 11)) ne 0; // not congruent

// ----------------
// (iii): The case when (N,r) = (3,1)
// J = J' case
E1 := EllipticCurve("14a5"); J1 := jInvariant(E1)/1728;
E2 := EllipticCurve("126a6"); J2 := jInvariant(E2)/1728;

Pa,Pb := CongPolys(3,1);
alpha := Roots(Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), x, 0]))[1][1];
assert #Roots(Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), alpha, x])) gt 0;

// Use the Hesse polynomials to show that E1 and E2 are not congruent
_,c4,c6 := HessePolynomials(3,1,cInvariants(E1));
f := c4^3 - J2*(c4^3 - c6^2);
pts := ProjectiveRoots(f);
cong_ell := [EllipticCurve([-27*Evaluate(c4,p), -54*Evaluate(c6,p)]) : p in pts];
assert not exists{E : E in cong_ell | IsIsomorphic(E2, E)};

// ----------------
// (iii): The case when (N,r) = (3,2)
E1 := EllipticCurve("245a1"); J1 := jInvariant(E1)/1728;
E2 := EllipticCurve("245b1"); J2 := jInvariant(E2)/1728;
Pa,Pb := CongPolys(3,2);
alpha := Roots(Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), x, 0]))[1][1];
assert #Roots(Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), alpha, x])) gt 0;

// Use the Hesse polynomials to show that E1 and E2 are not congruent
_,c4,c6 := HessePolynomials(3,2,cInvariants(E1));
f := c4^3 - J2*(c4^3 - c6^2);
pts := ProjectiveRoots(f);
assert #pts eq 0;

// ----------------
// (iii): The case when (N,r) = (4,1)
E1 := EllipticCurve("14a1"); J1 := jInvariant(E1)/1728;
E2 := EllipticCurve("112c3"); J2 := jInvariant(E2)/1728;
Pa,Pb,Pg := CongPolys(4,1);
assert exists(alpha){a[1] : a in Roots(Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), x, 0, 0])) | a[1] gt 0};
beta := Roots(Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), alpha, x, 0]))[1][1];
assert #Roots(Evaluate(Pg, [J1*J2, (J1-1)*(J2-1), alpha, beta, x])) gt 0;

// Use the Hesse polynomials to show that E1 and E2 are not congruent
_,c4,c6 := HessePolynomials(4,1,cInvariants(E1));
f := c4^3 - J2*(c4^3 - c6^2);
pts := ProjectiveRoots(f);
cong_ell := [EllipticCurve([-27*Evaluate(c4,p), -54*Evaluate(c6,p)]) : p in pts];
assert not exists{E : E in cong_ell | IsIsomorphic(E2, E)};

// ----------------
// (iii): The case when (N,r) = (4,1)
E3 := EllipticCurve("784j3");

// Follows from the (4,1) case
assert IsIsomorphic(E3, QuadraticTwist(E2, Discriminant(E2)));
