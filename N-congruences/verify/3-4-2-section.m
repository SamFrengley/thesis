/*
    Computational verifications for Z(15,11)
*/
Attach("../ZNr-equations.m");
QQ := Rationals();

// ---------------------------
// W(5, 1)
KK<a,b> := FunctionField(QQ, 2);

A := (16*a^2 + 208*a*b^2 - 312*a*b - 8*a + b^4 - 228*b^3 + 494*b^2 + 228*b + 1);

B := (64*a^3 + 4488*a^2*b^2 + 2448*a^2*b - 48*a^2 + 552*a*b^4 - 14076*a*b^3 - 
4692*a*b^2 - 2484*a*b + 12*a - b^6 - 522*b^5 + 10005*b^4 + 10005*b^2 + 522*b - 1);

JJ := A^(3)/((2)^(12)*(3)^(6)*(b)^(2)*(a)^(5));
JJ_1 := B^(2)/((2)^(12)*(3)^(6)*(b)^(2)*(a)^(5));
                  
c6_prod := -(64*a^3 + 4488*a^2*b^2 + 2448*a^2*b - 48*a^2 + 552*a*b^4 - 14076*a*b^3 - 
4692*a*b^2 - 2484*a*b + 12*a - b^6 - 522*b^5 + 10005*b^4 + 10005*b^2 + 522*b - 1);

//---------------------------
// W(5, 1)^{cbrt}
P<alpha_prime> := PolynomialRing(KK);

alpha := A/(2^4*3^2*a^2*b) * alpha_prime;

f := alpha^3 - JJ;
f := f * (2^12*3^6*b^3*a^6)/A^3; // rescale
assert f eq alpha_prime^3 - a*b;

//parametrise the surface
A3<a,b,alpha_prime> := AffineSpace(QQ, 3);
A2<c,d> := AffineSpace(QQ, 2);

W5_1_cbrt := Surface(A3, eval Sprint(f));
phi := map<A2 -> W5_1_cbrt | [ 1/(c^2*d), d/c, 1/c ]>;

assert IsInvertible(phi);

//---------------------------
//W+(15, 3)
KK<c,d> := FunctionField(QQ, 2);
P<eta_prime,theta_prime> := PolynomialRing(KK,2);

a,b,_ := Explode([ 1/(c^2*d), d/c, 1/c ]);
JJ := Evaluate(JJ, [a,b]); JJ_1 := Evaluate(JJ_1, [a,b]);
c6_prod := Evaluate(c6_prod, [a,b]);
_, alpha := IsPower(JJ,3);
B := Evaluate(B, [a,b]);

g1 := -1*(d)^(1)*(c^3*d - 119*c^2*d^2 + 12*c^2*d + 187*c*d^3 - 132*c*d^2 - 4*c - 
17*d^4 - 12*d^3 + 32*d - 12)^(1);
g2 := 1*(2)^(1)*(59*c^3*d^2 + 714*c^2*d^3 - 265*c^2*d^2 + 7*c^2*d + 323*c*d^4 + 
990*c*d^3 - 462*c*d^2 - 68*c*d + 12*d^5 - 85*d^4 - 77*d^3 - 576*d^2 + 160*d - 
28)^(1);
g3 := 1*(2)^(4)*(3)^(1)*(d)^(1)*(c^4*d^2 - 22*c^3*d^3 + 119*c^2*d^4 - 3*c^2*d + 
22*c*d^5 + 8*c*d^2 + d^6 - 47*d^3 + 1)^(1);
g4 := 1*(d)^(1)*(c^5*d^2 + 235*c^4*d^3 + 24*c^4*d^2 + 2090*c^3*d^4 - 1848*c^3*d^3 - 
8*c^3*d + 3230*c^2*d^5 + 17136*c^2*d^4 - 368*c^2*d^2 - 120*c^2*d - 1595*c*d^6 + 
4488*c*d^5 - 1976*c*d^3 + 1200*c*d^2 + 144*c*d + 16*c + 7*d^7 + 264*d^6 + 
1456*d^4 - 9720*d^3 + 1008*d^2 + 112*d + 96)^(1);
g5 := 1*(2)^(1)*(25*c*d - 35*d - 7)^(1)*(c^4*d^2 + 228*c^3*d^3 - 72*c^3*d^2 + 
494*c^2*d^4 + 768*c^2*d^3 - 8*c^2*d - 228*c*d^5 + 336*c*d^4 - 312*c*d^2 + 
168*c*d + d^6 + 24*d^5 + 208*d^3 - 744*d^2 + 144*d + 16)^(1);
g6 := 1*(d)^(1)*(c + 7*d)^(1);
g7 := 1*(2)^(1)*(25*c*d - 35*d - 7)^(1);

eta := (g1*theta_prime + g2)/(g6*theta_prime + g7);
theta := (g3*eta_prime + g4*theta_prime + g5)/(g6*theta_prime + g7);
beta := -B*c^6*d/((2)^(8)*(3)^(4)) * eta;

Q1 := eta_prime^2 + 4*eta_prime*theta_prime + (-3*c*d - d^2 + 4*d)*theta_prime^2 + (64*c*d + 28*c - 2*d^2 + 36*d - 
    64)*theta_prime - 49*c^2 + 214*c*d - 364*c - d^2 + 32*d - 256;
Q2 := (-2*c*d - 14*d^2)*eta_prime*theta_prime + (-100*c*d + 140*d + 28)*eta_prime + (-20*c*d^2 + 10*d^3 - 60*d^2 + 
    6*d)*theta_prime^2 + (-14*c^2*d - 236*c*d^2 + 32*c*d + 34*d^3 - 176*d^2 + 880*d)*theta_prime + 476*c^2*d + 
    764*c*d^2 + 2580*c*d + 196*c + 24*d^3 - 340*d^2 + 692*d - 3976;

f := beta^4 - 6*(alpha + 1)*JJ_1*beta^2 - 8*JJ_1^2*beta - 3*(alpha - 1)^2*JJ_1^2;
assert Numerator(f) in ideal<P | [Q1,Q2]>;
assert Numerator(eta^2 - theta) in ideal<P | [Q1,Q2]>;

// for later use need to store beta
P<eta_prime, theta_prime> := FunctionField(KK,2);
beta := Evaluate(beta, [eta_prime, theta_prime]);

//parametrise the surface
A4<c,d,eta_prime,theta_prime> := AffineSpace(QQ, 4);
A2<s,t> := AffineSpace(QQ, 2);

W15_11_plus := Surface(A4, [eval Sprint(Q1), eval Sprint(Q2)]);
phi := map<A2 -> W15_11_plus | [
-1*(s^3 - 2*s^2 - s*t^2 - s*t + s - t^3 + t^2 + t)*(s^3 - s^2 - s*t^2 - s*t -
t^2)/((s - t - 1)^(1)*(s*t - s + 1)*(s^2 - s*t^2 + s*t - s + t^3 - 2*t^2 - 2*t)),
-1*(s - t - 1)^(2)*(s - t + 3)/((s^2 - s*t^2 + s*t - s + t^3 - 2*t^2 - 2*t)^(1)),
-1*(7*s^8 + 7*s^7*t - 21*s^7 - 154*s^6*t^2 + 194*s^6*t - 68*s^6 + 149*s^5*t^3 +
20*s^5*t^2 - 672*s^5*t + 473*s^5 + 79*s^4*t^4 - 86*s^4*t^3 + 310*s^4*t^2 +
710*s^4*t - 981*s^4 - 108*s^3*t^5 + 44*s^3*t^4 + 323*s^3*t^3 - 103*s^3*t^2 -
398*s^3*t + 953*s^3 + 20*s^2*t^6 - 160*s^2*t^5 - 65*s^2*t^4 - 361*s^2*t^3 +
164*s^2*t^2 + 429*s^2*t - 438*s^2 + 72*s*t^6 - 213*s*t^5 - 241*s*t^4 + 125*s*t^3
- 267*s*t^2 - 357*s*t + 75*s + 7*t^6 + 44*t^5 - 161*t^4 - 206*t^3 + 30*t^2 +
87*t)/((s - t - 1)^(1)*(s*t - s + 1)*(s^2 - 4*s*t + 5*s + t - 6)*(s^2 - s*t^2 +
s*t - s + t^3 - 2*t^2 - 2*t)),
-1*(2)*(4*s^3 - 6*s^2*t + 2*s^2 - 19*s*t^2 + 14*s*t - 9*s - 4*t^2 - 15*t +
3)/((s - t - 1)^(1)*(s^2 - 4*s*t + 5*s + t - 6))
]>;

//assert IsInvertible(phi);

/* TODO
Should give explicit inverse in the nicest presentation you can
*/


//---------------------------
//W(15, 11)
KK<s,t> := FunctionField(QQ, 2);

c,d,eta_prime,theta_prime := Explode([
-1*(s^3 - 2*s^2 - s*t^2 - s*t + s - t^3 + t^2 + t)*(s^3 - s^2 - s*t^2 - s*t -
t^2)/((s - t - 1)^(1)*(s*t - s + 1)*(s^2 - s*t^2 + s*t - s + t^3 - 2*t^2 - 2*t)),
-1*(s - t - 1)^(2)*(s - t + 3)/((s^2 - s*t^2 + s*t - s + t^3 - 2*t^2 - 2*t)^(1)),
-1*(7*s^8 + 7*s^7*t - 21*s^7 - 154*s^6*t^2 + 194*s^6*t - 68*s^6 + 149*s^5*t^3 +
20*s^5*t^2 - 672*s^5*t + 473*s^5 + 79*s^4*t^4 - 86*s^4*t^3 + 310*s^4*t^2 +
710*s^4*t - 981*s^4 - 108*s^3*t^5 + 44*s^3*t^4 + 323*s^3*t^3 - 103*s^3*t^2 -
398*s^3*t + 953*s^3 + 20*s^2*t^6 - 160*s^2*t^5 - 65*s^2*t^4 - 361*s^2*t^3 +
164*s^2*t^2 + 429*s^2*t - 438*s^2 + 72*s*t^6 - 213*s*t^5 - 241*s*t^4 + 125*s*t^3
- 267*s*t^2 - 357*s*t + 75*s + 7*t^6 + 44*t^5 - 161*t^4 - 206*t^3 + 30*t^2 +
87*t)/((s - t - 1)^(1)*(s*t - s + 1)*(s^2 - 4*s*t + 5*s + t - 6)*(s^2 - s*t^2 +
s*t - s + t^3 - 2*t^2 - 2*t)),
-1*(2)*(4*s^3 - 6*s^2*t + 2*s^2 - 19*s*t^2 + 14*s*t - 9*s - 4*t^2 - 15*t +
3)/((s - t - 1)^(1)*(s^2 - 4*s*t + 5*s + t - 6))
]);

c6_prod := Evaluate(c6_prod, [c,d]);
beta := eval Sprint(beta);

assert IsSquare(3*c6_prod*beta / ((s - t - 1)*(s - t + 3)) );

//parametrise the surface
A3<s,t,w> := AffineSpace(QQ, 3);
A2<u,v> := AffineSpace(QQ, 2);

W15_11 := Surface(A3, w^2 - (s - t - 1)*(s - t + 3));
phi := map<A2 -> W15_11 | [
1*(u^2*v + u^2 + 2*u*v^2 + 2*u*v + 2*u + v^3 + v^2 + 2*v + 1)/((v + 1)^(1)*(u + 
v)*(u + v + 1)),
-1*(u*v - u + v^2)/((v + 1)^(1)*(u + v)*(u + v + 1)),
(2*u + 2*v + 1)/((u + v)*(u + v + 1))
]>;
assert IsInvertible(phi);

//---------------------------
//Z(15, 11)
KK<u,v> := FunctionField(QQ,2);
s,t,_ := Explode([
1*(u^2*v + u^2 + 2*u*v^2 + 2*u*v + 2*u + v^3 + v^2 + 2*v + 1)/((v + 1)^(1)*(u + 
v)*(u + v + 1)),
-1*(u*v - u + v^2)/((v + 1)^(1)*(u + v)*(u + v + 1)),
(2*u + 2*v + 1)/((u + v)*(u + v + 1))
]);

// Evaluating is faster if we do it in one step.
c,d := Explode([Evaluate(p, [s,t]) : p in [c,d]]);
JJ := Evaluate(JJ, [c,d]);
JJ_1 := Evaluate(JJ_1, [c,d]);

////////////////////
// Equation
f := -ZNrEquations(15,11 : vars:=[u,v,0])[1];

J_diff := (-JJ_1 + JJ + 1)^2 - 4*JJ;
assert IsSquare(f*J_diff);

////////////////////
// j-invs
jj, jj_1728 := WNrModuli(15,11 : vars:=[u,v]);
assert 1728^2*JJ eq jj;
assert 1728^2*JJ_1 eq jj_1728;

