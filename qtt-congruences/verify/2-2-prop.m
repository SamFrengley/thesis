/*
    The computations required to prove Proposition 4.2.2.
*/

load "../functions/XEeqns.m";

// -------------
// Define curves with CM j \neq 0, 1728 from Silverman's ``Advanced Topics
// in the Arithmetic of Elliptic Curves'' Appendix A, Section 3.
// -------------

E12:=EllipticCurve([0,0,0,-15,22]);
E27:=EllipticCurve([0,0,1,-30,63]);
E16:=EllipticCurve([0,0,0,-11,14]);
E7:=EllipticCurve([1,-1,0,-2,-1]);
E28:=EllipticCurve([0,0,0,-595,5586]);
E8:=EllipticCurve([0,4,0,2,0]);
E11:=EllipticCurve([0,-1,1,-7,10]);
E19:=EllipticCurve([0,0,1,-38,90]);
E43:=EllipticCurve([0,0,1,-860,9707]);
E67:=EllipticCurve([0,0,1,-7370,243528]);
E163:=EllipticCurve([0,0,1,-2174420,1234136692]);

CM:=[E12,E27,E16,E7,E28,E8,E11,E19,E43,E67,E163];

//---------------
// Claim D=8,11 traces of frobenius rule out 9 congruences.
// --------------

ell := 17;
E8_3 := QuadraticTwist(E8, -3);
E8_3D := QuadraticTwist(E8, 3*8);
E11_3 := QuadraticTwist(E11, -3);
E11_3D := QuadraticTwist(E11, 3*11);

// all have good reduction 
assert (Conductor(E8)*Conductor(E8_3)*Conductor(E8_3D)*Conductor(E11)*Conductor(E11_3)*Conductor(E11_3D) mod ell) ne 0;

assert (TrFrobDiff(E8, E8_3, ell) mod 9) ne 0;
assert (TrFrobDiff(E8, E8_3D, ell) mod 9) ne 0;
assert (TrFrobDiff(E8, E11_3, ell) mod 9) ne 0;
assert (TrFrobDiff(E8, E11_3D, ell) mod 9) ne 0;

//---------------
// Claim j=0 p-congruences are not p^2-congruences via traces of frobenius
// --------------

ell := 13;
E := EllipticCurve([0,10]);
E_5 := QuadraticTwist(E, 5);
E_15 := QuadraticTwist(E, -15);

assert (TrFrobDiff(E, E_5, ell) mod 25) ne 0;
assert (TrFrobDiff(E, E_15, ell) mod 25) ne 0;


E := EllipticCurve([0,98]);
E_7 := QuadraticTwist(E, -7);
E_21 := QuadraticTwist(E, 21);

assert (TrFrobDiff(E, E_7, ell) mod 49) ne 0;
assert (TrFrobDiff(E, E_21, ell) mod 49) ne 0;

//---------------
// Proof that the 4-congruences in Table 4.1 are the only 4-congruences 
// between elliptic curves with CM and j \neq 0, 1728.
// --------------

// Define Table 4.1 4-congruences with j \neq 0, 1728.
table41 := {<E7, QuadraticTwist(E7, -1),1>,
            <E7, QuadraticTwist(E7, 7),3>,
            <E8, QuadraticTwist(E8, 2),3>,
            <E12, QuadraticTwist(E12, -1),1>,
            <E12, QuadraticTwist(E12, 3),3>,
            <E16, QuadraticTwist(E16, 2),3>,
            <E28, QuadraticTwist(E28, -1),1>,
            <E28, QuadraticTwist(E28, 7),3> };

table41 := {<info[1], MinimalModel(info[2]), info[3]> : info in table41};

// Determine nontrivially (4,r)-congruent curves to the CM curves 
our_table41 := {};
for E in CM do
    cc := cInvariants(E);
    ab := [-27*cc[1], -54*cc[2]];
    EE1 := NonTrivially4CongruentQT(ab, 1);
    EE3 := NonTrivially4CongruentQT(ab, 3);
    our_table41 := our_table41 join {<E, MinimalModel(E1[1]), 1> : E1 in EE1};
    our_table41 := our_table41 join {<E, MinimalModel(E3[1]), 3> : E3 in EE3};
end for;

// Check the tables are equal
assert our_table41 eq table41;

//---------------
// Proof that the 4-congruences in Table 4.1 where j \neq 0, 1728 are not 8-congruences
// --------------

for cong in table41 do
    E := cong[1];
    cc := cInvariants(E);
    ab := [-27*cc[1], -54*cc[2]];

    assert #NonTrivially8CongruentQT(ab, 1) eq 0;
    assert #NonTrivially8CongruentQT(ab, 3) eq 0;
    assert #NonTrivially8CongruentQT(ab, 5) eq 0;
    assert #NonTrivially8CongruentQT(ab, 7) eq 0;
end for;

//---------------
// 4-congruences with j=0
// --------------
Qa<a> := FunctionField(Rationals());
XEr := ProjectiveSpace(Qa, 1);
cc := [0, -a/54];

_, c4, c6 := HessePolynomials(4, 1, cc : Variables:=[XEr.1,XEr.2]);
fcts := Factorisation(c4);
assert fcts[1][1] eq XEr.2;
assert fcts[2][1] eq XEr.1;
assert 27*fcts[3][1] eq 27*XEr.1^3 - 8*a*XEr.2^3;
assert 27*fcts[4][1] eq 27*XEr.1^3 + 1*a*XEr.2^3;

// since a cubefree >0 the last two only have roots when a=1
// for the first two just show it doesn't work unless a=1
aa2 := [-27*Evaluate(c4, [1,0]), -54*Evaluate(c6, [1,0])];
aa3 := [-27*Evaluate(c4, [0,1]), -54*Evaluate(c6, [0,1])];

assert IsIsomorphic(EllipticCurve([0,a]), EllipticCurve(aa2));
// the trivial quadratic twist. The discriminant is -3 up to a 
// square so this is trivial for (4,3) too.

assert aa3[2]/a eq -(2^6/3^(12))*a^4;
// this is not a cube unless a=1 so is not a quadratic twist

// so it all reduces to the case when a=1. In this case

j0_41 := NonTrivially4CongruentQT([0,1], 1);
j0_43 := NonTrivially4CongruentQT([0,1], 3);

for E2 in j0_41 do 
    assert IsIsomorphic(E2[1], EllipticCurve([0,-1]));
end for;
for E2 in j0_43 do
    assert IsIsomorphic(E2[1], EllipticCurve([0,27]));
end for;

//---------------
// Claim j=0 4-congruences are not 8-congruences
//---------------

E := EllipticCurve([0,1]);
assert #NonTrivially8CongruentQT([0,1], 1) eq 0;
assert #NonTrivially8CongruentQT([0,1], 3) eq 0;
assert #NonTrivially8CongruentQT([0,1], 5) eq 0;
assert #NonTrivially8CongruentQT([0,1], 7) eq 0;


//---------------
// 4 and 8-congruences with j=1728
// --------------
Qa<a> := FunctionField(Rationals());
P<x> := PolynomialRing(Qa);
XEr := ProjectiveSpace(Qa, 1);
cc := [-a/27, 0];

_, c4, c6 := HessePolynomials(4, 1, cc : Variables:=[XEr.1,XEr.2]);
fcts := Factorisation(c6);
assert fcts[1][1] eq XEr.2;
assert fcts[2][1] eq XEr.1;
assert 9*fcts[3][1] eq 9*XEr.1^2 + a*XEr.2^2;
assert 9*fcts[4][1] eq 9*XEr.1^4 - 2*a*XEr.1^2*XEr.2^2 + a^2*XEr.2^4;
assert 729*fcts[5][1] eq 729*XEr.1^4 - 18*a*XEr.1^2*XEr.2^2 + a^2*XEr.2^4;

// since a squarefree the third only has roots when a=1
// the last two never have points since the discriminant of the 
// quadratic is negative.
aa2 := [-27*Evaluate(c4, [1,0]), -54*Evaluate(c6, [1,0])];
aa3 := [-27*Evaluate(c4, [0,1]), -54*Evaluate(c6, [0,1])];

assert IsIsomorphic(EllipticCurve([a,0]), EllipticCurve(aa2));
assert IsIsomorphic(EllipticCurve([a,0]), EllipticCurve(aa3));
// the trivial quadratic twist. 
// The (4,3) congruent one is the quadratic twist by the
// discriminant which is up to a square equal to a.

// The case when a=1.
j1728_41 := NonTrivially4CongruentQT([1,0], 1);
j1728_43 := NonTrivially4CongruentQT([1,0], 3);
assert #j1728_41 eq 0;
assert #j1728_43 eq 0;

// it is elementary to check by hand that the points [0,1] and [1,0] only 
// lift to X_E^3(8) and X_E^7(8) when a = +2 and -2 respectively.

//---------------
// Claim j=1728 8-congruences are not 16-congruences via trace of frobenius
//---------------

ell := 29;
E := EllipticCurve([2,0]);
Ed := QuadraticTwist(E, 2);

assert (TrFrobDiff(E, Ed, ell) mod 16) ne 0;

E := EllipticCurve([-2,0]);
Ed := QuadraticTwist(E, 2);

assert (TrFrobDiff(E, Ed, ell) mod 16) ne 0;
