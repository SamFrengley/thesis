/*
    Computational verifications for Z(14,1)
*/
load "../ZNr-equations.m";
QQ := Rationals();

//---------------------------
//W(7, 1)
KK<a,b> := FunctionField(QQ, 2);

A := 2^4*(a^4*b^2 - 2*a^4*b + a^4 + 44*a^3*b^3 + 64*a^3*b^2 - 52*a^3*b - 
    56*a^3 + 21*a^2*b^4 + 23*a^2*b^3 + 11*a^2*b^2 - 27*a^2*b - 12*a^2 + 14*a*b^5 + 
    7*a*b^4 - 39*a*b^3 - 4*a*b^2 + 25*a*b + 13*a + 4*b^6 - 2*b^5 - b^4 - 9*b^3 - 
    3*b^2 + 11*b + 4);

B := (64*a^6*b^3 - 192*a^6*b^2 + 192*a^6*b - 64*a^6 - 9600*a^5*b^4 + 1920*a^5*b^3 +
    16512*a^5*b^2 - 384*a^5*b - 8448*a^5 + 7008*a^4*b^5 - 46272*a^4*b^4 - 
    79872*a^4*b^3 + 12864*a^4*b^2 + 76704*a^4*b + 29568*a^4 + 6496*a^3*b^6 - 
    29856*a^3*b^5 - 155232*a^3*b^4 - 107616*a^3*b^3 + 104256*a^3*b^2 + 140544*a^3*b 
    + 37312*a^3 - 10056*a^2*b^7 - 2136*a^2*b^6 - 38064*a^2*b^5 - 134736*a^2*b^4 - 
    45768*a^2*b^3 + 113064*a^2*b^2 + 94656*a^2*b + 16896*a^2 - 3444*a*b^8 - 
    6396*a*b^7 + 1116*a*b^6 - 19932*a*b^5 - 44508*a*b^4 - 2964*a*b^3 + 39540*a*b^2 +
    29292*a*b + 4224*a - 539*b^9 - 21*b^8 - 1584*b^7 - 268*b^6 - 1974*b^5 - 5526*b^4
    + 608*b^3 + 4764*b^2 + 3489*b + 539);

JJ := A^(3)/((3)^(6)*(b - 1)*(b + 1)^(3)*(8*a*b + b^2 + 4*b - 1)^(7));

JJ_1 := B^(2)/((3)^(6)*(b - 1)*(b + 1)^(3)*(8*a*b + b^2 + 4*b - 1)^(7));

//---------------------------
//W(7, 1)^{sqrt}
P<alpha_prime> := PolynomialRing(KK);

alpha := B/(3^3*(b+1)^2*(8*a*b + b^2 + 4*b - 1)^4) * alpha_prime;

f := alpha^2 - JJ_1;
f := f * (3^6*(b+1)^4*(b-1)*(8*a*b + b^2 + 4*b - 1)^8)/B^2; // rescale
assert f eq (b - 1)*alpha_prime^2 - (b+1)*(8*a*b + b^2 + 4*b - 1);

//parametrise the surface
A3<a,b,alpha_prime> := AffineSpace(QQ, 3);
A2<c,d> := AffineSpace(QQ, 2);

W7_1_sqrt := Surface(A3, eval Sprint(f));
phi := map<A2 -> W7_1_sqrt | [-1*(c^2*d - c^2 - d^3 + 3*d^2 + 5*d + 1)/((2)^(3)*(d)*(d + 1)), 1/d, c/d]>;
assert IsInvertible(phi);

//---------------------------
//W(14, 1)
KK<c,d> := FunctionField(QQ, 2);
P<beta_prime> := PolynomialRing(KK);

a,b,_ := Explode([-1*(c^2*d - c^2 - d^3 + 3*d^2 + 5*d + 1)/((2)^(3)*(d)*(d + 1)), 1/d, c/d]);
JJ := Evaluate(JJ, [a,b]); JJ_1 := Evaluate(JJ_1, [a,b]);
_, alpha := IsSquare(JJ_1);

h1 := -1*((d^3 - 3*d^2 + 3*d - 1)*c^7 + (-d^4 - 36*d^3 + 34*d^2 + 44*d - 41)*c^6 + 
    (-3*d^5 - 225*d^4 - 350*d^3 + 102*d^2 + 353*d + 123)*c^5 + (3*d^6 + 342*d^5 + 
    341*d^4 - 276*d^3 + 181*d^2 + 958*d + 499)*c^4 + (3*d^7 + 459*d^6 + 895*d^5 + 
    1671*d^4 + 3497*d^3 + 1985*d^2 - 1579*d - 1299)*c^3 + (-3*d^8 - 576*d^7 - 
    988*d^6 - 3072*d^5 - 9866*d^4 - 10048*d^3 + 948*d^2 + 6016*d + 2229)*c^2 + (-d^9
    - 231*d^8 - 548*d^7 - 2852*d^6 - 8878*d^5 - 9962*d^4 - 436*d^3 + 7388*d^2 + 
    5511*d + 1305)*c + d^10 + 270*d^9 + 613*d^8 + 4384*d^7 + 16554*d^6 + 24684*d^5 +
    11770*d^4 - 10224*d^3 - 18187*d^2 - 10922*d - 2559);

h2 := -1*((d^3 - 3*d^2 + 3*d - 1)*c^7 + (-17*d^4 - 6*d^3 + 40*d^2 + 6*d - 23)*c^6 + 
    (-3*d^5 - 189*d^4 - 334*d^3 + 14*d^2 + 337*d + 175)*c^5 + (51*d^6 - 360*d^5 - 
    805*d^4 + 400*d^3 + 1649*d^2 + 984*d + 129)*c^4 + (3*d^7 + 387*d^6 + 1391*d^5 + 
    2159*d^4 + 1929*d^3 + 905*d^2 + 5*d - 123)*c^3 + (-51*d^8 + 738*d^7 + 1718*d^6 -
    4118*d^5 - 14480*d^4 - 11482*d^3 + 2410*d^2 + 6670*d + 2211)*c^2 + (-d^9 - 
    195*d^8 - 1060*d^7 - 3108*d^6 - 4174*d^5 - 930*d^4 + 2188*d^3 + 188*d^2 - 1817*d
    - 819)*c + 17*d^10 - 372*d^9 - 953*d^8 + 4648*d^7 + 19442*d^6 + 22240*d^5 - 
    394*d^4 - 17896*d^3 - 10419*d^2 + 84*d + 1011);

h3 := (d - 1)*c^3 + (-d^2 + 10*d + 7)*c^2 + (-d^3 + 19*d^2 + 29*d + 9)*c + d^4 - 
    28*d^3 - 66*d^2 - 36*d + 1;

h4 := (d - 1)*c^3 + (7*d^2 + 8*d + 1)*c^2 + (-d^3 + 7*d^2 + 21*d + 13)*c - 7*d^4 - 
    14*d^3 - 28*d^2 - 42*d - 21;

beta := (d)^(6)*Evaluate(A, [a,b])/(2^4*(3)^(3)*(d - 1)^(3)*(c)^(7))*(h1*beta_prime + h2)/(h3*beta_prime + h4);
f := Numerator(beta^3 - 3*JJ*beta - 2*JJ*(alpha + 1));
f := f * LCM([Denominator(c) : c in Coefficients(f)]) * 782757789696; // scale
f := f / GCD([Numerator(c) : c in Coefficients(f)]); // scale

//parametrise the surface
A3<c,d,beta_prime> := AffineSpace(QQ, 3);
A2<u,v> := AffineSpace(QQ, 2);

W14_1 := Surface(A3, eval Sprint(f));
phi := map<A2 -> W14_1 | [
    -2*(v)*(u - 1)*(u^2*v - u^2 - u*v^2 + u*v - u + v^2)/((u*v - u - v - 1)*(u^2 - u*v + u + v)), 
    -1*(u^3*v - u^3 + u^2*v^2 + u^2*v - 2*u^2 - 2*u*v^2 - u*v - u + v^2 - v)/((u*v -u - v - 1)^(1)*(u^2 - u*v + u + v)), 
    (-u*v + v)/(u*v - u - v - 1)]>;
assert IsInvertible(phi);

//---------------------------
//W(14, 1)
KK<u,v> := FunctionField(QQ, 2);

c,d,_ := Explode([
    -2*(v)*(u - 1)*(u^2*v - u^2 - u*v^2 + u*v - u + v^2)/((u*v - u - v - 1)*(u^2 - u*v + u + v)), 
    -1*(u^3*v - u^3 + u^2*v^2 + u^2*v - 2*u^2 - 2*u*v^2 - u*v - u + v^2 - v)/((u*v -u - v - 1)^(1)*(u^2 - u*v + u + v)), 
    (-u*v + v)/(u*v - u - v - 1)]
);

JJ := Evaluate(JJ, [c,d]); JJ_1 := Evaluate(JJ_1, [c,d]);

////////////////////
// Equation
f := -ZNrEquations(14,1 : vars:=[u,v,0])[1];

J_diff := (-JJ_1 + JJ + 1)^2 - 4*JJ;
assert IsSquare(f*J_diff);


////////////////////
// j-invs
jj, jj_1728 := WNrModuli(14,1 : vars:=[u,v]);
assert 1728^2*JJ eq jj;
assert 1728^2*JJ_1 eq jj_1728;
