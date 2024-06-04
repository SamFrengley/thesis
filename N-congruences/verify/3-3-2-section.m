/*
    Computational verifications for Z(14,-1)
*/
load "../ZNr-equations.m";
QQ := Rationals();

//---------------------------
//W(7, -1)
KK<a,b> := FunctionField(QQ, 2);

A := (144*a^4*b^6 - 48*a^4*b^5 - 80*a^4*b^4 + 320*a^4*b^3 + 200*a^4*b^2 - 8*a^4*b +
    a^4 + 384*a^3*b^5 + 824*a^3*b^4 + 528*a^3*b^3 - 136*a^3*b^2 - 220*a^3*b + 
    216*a^2*b^5 + 632*a^2*b^4 + 608*a^2*b^3 + 198*a^2*b^2 + 16*a*b^4 + 20*a*b^3 + 
    b^4);
B := (12960*a^6*b^8 + 31680*a^6*b^7 + 25768*a^6*b^6 + 5736*a^6*b^5 - 2880*a^6*b^4 -
    1040*a^6*b^3 + 540*a^6*b^2 + 12*a^6*b - a^6 - 2592*a^5*b^8 - 7056*a^5*b^7 - 
    25440*a^5*b^6 - 47184*a^5*b^5 - 32292*a^5*b^4 - 4080*a^5*b^3 + 1476*a^5*b^2 - 
    534*a^5*b + 6480*a^4*b^8 + 28728*a^4*b^7 + 47520*a^4*b^6 + 42804*a^4*b^5 + 
    34512*a^4*b^4 + 24468*a^4*b^3 + 7473*a^4*b^2 + 4608*a^3*b^7 + 16708*a^3*b^6 + 
    22272*a^3*b^5 + 12372*a^3*b^4 + 2180*a^3*b^3 + 540*a^2*b^7 + 1548*a^2*b^6 + 
    1440*a^2*b^5 + 417*a^2*b^4 - 24*a*b^6 - 30*a*b^5 - 
    b^6);

JJ := A^(3)/((2)^(12)*(3)^(6)*(b)^(7)*(b + 1)^(4)*(a)^(5)*(a*b + a - 1)^(7));

JJ_1 := B^(2)/((2)^(12)*(3)^(6)*(b)^(7)*(b + 1)^(4)*(a)^(5)*(a*b + a - 1)^(7));

//---------------------------
//W(7, -1)^{sqrt}
P<alpha_prime> := PolynomialRing(KK);

alpha := B/(2^6*3^3*a^3*b^3*(b+1)^2*(a*b + a - 1)^4) * alpha_prime;

f := alpha^2 - JJ_1;
f := f * ((2)^(12)*(3)^(6)*(b)^(7)*(b + 1)^(4)*(a)^(6)*(a*b + a - 1)^(8))/B^2; // rescale
assert f eq b*alpha_prime^2 - a*(a*b + a - 1);

//parametrise the surface
A3<a,b,alpha_prime> := AffineSpace(QQ, 3);
A2<c,d> := AffineSpace(QQ, 2);

W7_6_sqrt := Surface(A3, eval Sprint(f));
phi := map<A2 -> W7_6_sqrt | [(-d + 1)/(c^2 - d^2 + d), d -1, c/(c^2 - d^2 + d)]>;
assert IsInvertible(phi);

//---------------------------
//W(14, 13)
KK<c,d> := FunctionField(QQ, 2);
P<beta_prime> := PolynomialRing(KK);

a,b,_ := Explode([(-d + 1)/(c^2 - d^2 + d), d - 1, c/(c^2 - d^2 + d)]);
JJ := Evaluate(JJ, [a,b]); JJ_1 := Evaluate(JJ_1, [a,b]);
_, alpha := IsSquare(JJ_1);

h1 := (-c^6 + 18*c^5*d + 2*c^5 - 15*c^4*d^2 - c^4*d - 15*c^4 - 36*c^3*d^3 - 40*c^3*d^2
    + 48*c^3*d - 4*c^3 + 33*c^2*d^4 + 98*c^2*d^3 - 253*c^2*d^2 - 54*c^2*d + 
    33*c^2 + 18*c*d^5 + 110*c*d^4 - 78*c*d^3 - 198*c*d^2 + 166*c*d + 2*c - 
    17*d^6 - 169*d^5 + 310*d^4 + 395*d^3 - 790*d^2 + 271*d - 17);

h2 := (c^7 - c^6*d - 3*c^5*d^2 - 23*c^5*d - 5*c^5 + 3*c^4*d^3 - 19*c^4*d^2 - 7*c^4*d + 
    18*c^4 + 3*c^3*d^4 - 38*c^3*d^3 + 113*c^3*d^2 - 18*c^3*d + 7*c^3 - 3*c^2*d^5
    + 122*c^2*d^4 - 317*c^2*d^3 + 642*c^2*d^2 - 19*c^2*d - 36*c^2 - c*d^6 + 
    61*c*d^5 - 408*c*d^4 + 289*c*d^3 + 280*c*d^2 - 239*c*d - 3*c + d^7 - 103*d^6
    + 696*d^5 - 1525*d^4 + 740*d^3 + 347*d^2 - 141*d + 18);

h3 := c^2 + (6*d - 2)*c - 7*d^2 + 7*d - 7;

h4 := -c^3 + d*c^2 + (d^2 - 7*d + 3)*c - d^3 + d^2 + 9*d + 6;

beta := (c^2 - d^2 + d)^(4)*Evaluate(A, [a,b])/((2)^(6)*(3)^(3)*(d)^(2)*(d - 1)^(4)*(c)^(7)) *  (h1 *beta_prime + h2) / (h3*beta_prime + h4);

f := Numerator(beta^3 - 3*JJ*beta - 2*JJ*(alpha + 1));
f := f * LCM([Denominator(c) : c in Coefficients(f)]) * 2985984; // scale
f := -f / GCD([Numerator(c) : c in Coefficients(f)]); // scale

assert f eq beta_prime^3 + (-c - d - 4)*beta_prime^2 + (2*c + 3*d + 5)*beta_prime + c*d - c  + d^2 - 3*d - 2;

//parametrise the surface
A3<c,d,beta_prime> := AffineSpace(QQ, 3);
A2<u,v> := AffineSpace(QQ, 2);

W14_13 := Surface(A3, eval Sprint(f));
phi := map<A2 -> W14_13 | [
    -(u^2 + 2*u*v^2 + 2*u - 2*v^3 + v^2 + 4*v + 1)/((2)^(2)*(v - 1)*(v + 1)), 
    (u - v - 1)*(u - v + 3)/((2)^(2)*(v - 1)*(v + 1)), 
    -1*(u - 3*v + 5)/((2)^(1)*(v - 1))]>;
assert IsInvertible(phi);


//---------------------------
//Z(14, 13)
KK<u,v> := FunctionField(QQ, 2);

c,d,_ := Explode([
    (u^2 + 2*u*v^2 + 2*u - 2*v^3 + v^2 + 4*v + 1)/((2)^(2)*(v - 1)*(v + 1)), 
    (u - v - 1)*(u - v + 3)/((2)^(2)*(v - 1)*(v + 1)), 
    -1*(u - 3*v + 5)/((2)^(1)*(v - 1))]
);

JJ := Evaluate(JJ, [c,d]); JJ_1 := Evaluate(JJ_1, [c,d]);

////////////////////
// Equation
f := -ZNrEquations(14,13 : vars:=[u,v,0])[1];

J_diff := (-JJ_1 + JJ + 1)^2 - 4*JJ;
assert IsSquare(f*J_diff);

////////////////////
// j-invs
jj, jj_1728 := WNrModuli(14,13 : vars:=[u,v]);
assert 1728^2*JJ eq jj;
assert 1728^2*JJ_1 eq jj_1728;
