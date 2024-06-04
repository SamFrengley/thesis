/*
    Computational verifications for Section 5.2.
*/
load "../Z12-r/Z12equations.m";
load "../Z2_3_4/polys.m";   

QQ := Rationals();

Pa_3, Pb_3 := CongPolys(3, 1);
Pa_4, Pb_4, Pg_4 := CongPolys(4, 1);

//---------------------------
//W(3, 1)
W31<a,b> := AffineSpace(QQ, 2);
_<x> := PolynomialRing(FieldOfFractions(CoordinateRing(W31)));

JJ := (2*a*b + 4*a - b^2 + 4*b)^(3)/((2*b + 1)^(3)*(a*b + 2*a - b^2 + 2*b)^(2));
JJ_1 := (b + 2)^(2)*(a^2*b + 2*a^2 - 2*a*b^2 + a*b + a + b^3 - 2*b^2 + b)^(2)/((2*b + 1)^(3)*(a*b + 2*a - b^2 + 2*b)^(2));

A4<alpha_3_var, beta_3_var, JJ_var, JJ_1_var> := AffineSpace(QQ, 4);

calW31 := Surface(A4, [ Evaluate(Pa_3, [JJ_var, JJ_1_var, alpha_3_var, beta_3_var]), 
                        Evaluate(Pb_3, [JJ_var, JJ_1_var, alpha_3_var, beta_3_var]) ]);


alpha_3 := Roots(Evaluate(Pa_3, [JJ, JJ_1, x, 0]))[1][1];
beta_3 := Roots(Evaluate(Pb_3, [JJ, JJ_1, alpha_3, x]))[1][1];

phi := map<W31 -> calW31 | [alpha_3, beta_3, JJ, JJ_1]>;
assert IsInvertible(phi);

//--------------------------
//W(3, 1)^{sqrt(DD')}
A3<a,b,alpha_4> := AffineSpace(QQ, 3); 
S1 := Surface(A3, Numerator(Evaluate(Pa_4, [Evaluate(JJ, [a,b]), Evaluate(JJ_1, [a,b]), alpha_4, 0, 0])));

A3<a,b,alpha_4_prime> := AffineSpace(QQ, 3); 
S2 := Surface(A3, alpha_4_prime^2 - (2*b + 1));

alpha_4 := -(b + 2)*(a^2*b - 2*a*b^2 + b^3 + 2*a^2 + a*b - 2*b^2 + a + b)/((2*b + 1)^2*(-a*b + b^2 - 2*a - 2*b))*alpha_4_prime;
phi := map<S2 -> S1 | [a,b,alpha_4]>;
assert IsInvertible(phi);

W31_sqrtDD<c,d> := AffineSpace(QQ, 2);
a := c;
b := (d^2 - 1)/2;
alpha_4_prime := d;

phi := map<W31_sqrtDD -> S2 | [a,b,alpha_4_prime]>;
assert IsInvertible(phi);

JJ := Evaluate(JJ, [a, b]);
JJ_1 := Evaluate(JJ_1, [a, b]);
alpha_3 := Evaluate(alpha_3, [a, b]);
beta_3 := Evaluate(beta_3, [a, b]);
alpha_4 := Evaluate(alpha_4, [a,b,alpha_4_prime]);

//--------------------------
//W(6, 1)
A3<c,d,beta_4> := AffineSpace(QQ, 3); 
S1 := Surface(A3, Numerator(Evaluate(Pb_4, [Evaluate(JJ, [c,d]), Evaluate(JJ_1, [c,d]), Evaluate(alpha_4, [c,d]), beta_4, 0])));

A3<c,d,beta_4_prime> := AffineSpace(QQ, 3); 
S2 := Surface(A3, beta_4_prime^3 + 6*beta_4_prime^2 - 3*(d + 1)*(d - 3)*beta_4_prime - (2*c*d^2 + 6*c - d^4 + 2*d^3 + 6*d^2 - 6*d - 9));

beta_4 := ((d^(4) - 4*c*d^(2) - 10*d^(2) - 12*c+ 9)*((d^(2) - 2*d - 3)*beta_4_prime - d^(4) + 2*c*d^(2) + 2*d^(3) + 6*d^(2) + 6*c - 6*d - 9 ))/(2*d^3*(d^(4) - 2*c*d^(2) - 6*d^(2) - 6*c + 5) )*1/beta_4_prime;
phi := map<S2 -> S1 | [c,d,beta_4]>;
assert IsInvertible(phi);

W61<p,q> := AffineSpace(QQ, 2);
c := (p^3*q + 6*p^2*q^2 + 9*p*q^3 + 12*p*q^2 - 12*p*q + 9*q^4 + 12*q^3 - 24*q^2 - 16*q + 16)/(2*q^2*(3*q^2 + 4));
d := 2/q;
beta_4_prime := p/q;

phi := map<W61 -> S2 | [c, d, beta_4_prime]>;
assert IsInvertible(phi);

JJ := Evaluate(JJ, [c, d]);
JJ_1 := Evaluate(JJ_1, [c, d]);
alpha_3 := Evaluate(alpha_3, [c, d]);
beta_3 := Evaluate(beta_3, [c, d]);
alpha_4 := Evaluate(alpha_4, [c, d]);
beta_4 := Evaluate(beta_4, [c, d, beta_4_prime]);

//--------------------------
//W(12, 1)^+
A3<p,q,gamma_4> := AffineSpace(QQ, 3); 
f1 :=  Evaluate(Pb_4, [Evaluate(JJ, [p,q]), Evaluate(JJ_1, [p,q]), Evaluate(alpha_4, [p,q]), Evaluate(beta_4, [p,q]), gamma_4]);

A3<p,q,gamma_4_prime> := AffineSpace(QQ, 3); 
f2 := (q + 1/4*p - 1)*gamma_4_prime^4 + (-9*q^2 - 9*q*p - 12*q - 3/2*p^2 + 12)*gamma_4_prime^2 + (-16*q - 16*p - 32)*gamma_4_prime - 3/4*p^3 - 3*p^2;

gamma_4 := (2*p^3*q + 12*p^2*q^2 + 18*p*q^3 + 24*p*q^2 - 24*p*q + 9*q^4 + 24*q^3 - 8*q^2 - 32*q + 16)^(2)/(64*(p + q + 2)^(3)*(p + 4*q - 4))*gamma_4_prime;
assert Numerator(Evaluate(f1, [p,q,gamma_4])) in Ideal(f2);

A3<p,q_prime,gamma_4_prime> := AffineSpace(QQ, 3); 
f3 := 3*q_prime^2 - 2*gamma_4_prime*q_prime + 8*q_prime + p + 4;

q := -(3*p*q_prime + 2*p*gamma_4_prime + 4*p - q_prime*gamma_4_prime^2 - 4*q_prime*gamma_4_prime + 8*q_prime - 4*gamma_4_prime + 16)/(6*(gamma_4_prime));
assert Numerator(Evaluate(f2, [p,q,gamma_4_prime])) in Ideal(f3);

// The above shows that the image of V(f3) -> A3 : (u',v',w) \mapsto (u, v, w) is V(f1).
// The map is invertible since you can (1) solve for v from (u, v', w') then
// (2) solve for w from (u, v w').

S := Surface(A3, f3);
W12_plus<s,t> := AffineSpace(QQ, 2);

p := -1*(2)^(2)*(s - 3)*(s + 1)*(3*s - 2*t - 1)^(2)/((3)^(1)*(s - 1)^(2)*(s - 2*t - 1)*(3*s + 2*t + 1));
q_prime := 1*(2)^(3)*(t - 1)*(s)/((3)^(1)*(s - 1)*(s - 2*t - 1));
gamma_4_prime := 1*(2)^(2)*(s)*(s - 3)*(3*s - 2*t - 1)/((s - 1)^(1)*(s - 2*t - 1)*(3*s + 2*t + 1));
q := Evaluate(q, [p, q_prime, gamma_4_prime]);

phi := map<W12_plus -> S | [p, q_prime, gamma_4_prime]>;
assert IsInvertible(phi);

JJ := Evaluate(JJ, [p, q]);
JJ_1 := Evaluate(JJ_1, [p, q]);
alpha_4 := Evaluate(alpha_4, [p, q]);
beta_4 := Evaluate(beta_4, [p, q]);
gamma_4 := Evaluate(gamma_4, [p, q, gamma_4_prime]);
alpha_3 := Evaluate(alpha_3, [p, q]);
beta_3 := Evaluate(beta_3, [p, q]);
JpJ := JJ - JJ_1 + 1;
delta_3 := -6*(2*beta_3^3 - (5*alpha_3 + 2)*beta_3^2 - 10*JJ_1*beta_3 + 3*JJ_1*(13*alpha_3 - 2 + 6*JpJ))/(beta_3^3 - 3*JJ_1*beta_3 - 2*JJ_1^2);
delta_4 := 3*(JJ - (alpha_4 + 1)^2);

//--------------------------
//W(12, 1)
A3<s,t,w> := AffineSpace(QQ, 3); 

f := 3*(t - 1)*(s)*(-s + 2*t + 1)*(3*s + 2*t + 1)*(s*t + 2*s - 2*t - 1);
assert IsSquare(Evaluate(alpha_4*delta_3*delta_4, [s,t])*f);
S := Surface(A3, w^2 - f);

W12_1<u,v> := AffineSpace(QQ, 2);

s := -(2*u^2 + u*v^2 + u - 4*v^2)/((u*v^2 + u + 4*v^2));
t := (2*u^2*v^2 - 4*u^2 + u*v^4 - 6*u*v^2 - 3*u - 8*v^4)/((2*u + v^2 + 3)*(u*v^2 + u + 4*v^2));
w := 36*(v)*(u + 2)*(u + v^2 + 1)*(u + 2*v^2)*(2*u^2 + u*v^2 + u - 4*v^2)^(2)/((2*u + v^2 + 3)^(2)*(u*v^2 + u + 4*v^2)^(3));

phi := map<W12_1 -> S | [s, t, w]>;
assert IsInvertible(phi);

JJ_W12_1 := Evaluate(JJ, [s, t]);
JJ_1_W12_1 := Evaluate(JJ_1, [s, t]);
JpJ_W12_1 := JJ_W12_1 - JJ_1_W12_1 + 1;

assert IsSquare((JpJ_W12_1^2 - 4*JJ_W12_1)*Z12Equations(1, u, v));

formulae_JJ, formulae_JJ_1 := W12Moduli(1, u, v);
assert 1728^2*JJ_W12_1 eq formulae_JJ;
assert 1728^2*JJ_1_W12_1 eq formulae_JJ_1;


//--------------------------
//W(12, 7)
A3<s,t,w> := AffineSpace(QQ, 3); 

f := -3*(t - 1)*(9*s^3 - 12*s^2*t - 15*s^2 - 4*s*t^2 - 4*s*t - s + 12*t^2 + 12*t + 3);
assert IsSquare(Evaluate(alpha_4*delta_3, [s,t])*f);
S := Surface(A3, w^2 - f);

W12_7<u,v> := AffineSpace(QQ, 2);

s := (u - 1)*(u^2 + 4*u - v^2 + 1)/((u^3 + u^2 - u*v^2 - 3*u - v^2 + 1));
t := (u^3 + 7*u^2 - u*v^2 + 9*u - v^2 + 1)/((u^3 + u^2 - u*v^2 - 3*u - v^2 + 1));
w := 36*(v)*(u)*(u + 2)*(u^2 + 4*u - v^2 + 1)/((u^3 + u^2 - u*v^2 - 3*u - v^2 + 1)^(2));

phi := map<W12_7 -> S | [s, t, w]>;
assert IsInvertible(phi);

JJ_W12_7 := Evaluate(JJ, [s, t]);
JJ_1_W12_7 := Evaluate(JJ_1, [s, t]);
JpJ_W12_7 := JJ_W12_7 - JJ_1_W12_7 + 1;

assert IsSquare((JpJ_W12_7^2 - 4*JJ_W12_7)*Z12Equations(7, u, v));

formulae_JJ, formulae_JJ_1 := W12Moduli(7, u, v);
assert 1728^2*JJ_W12_7 eq formulae_JJ;
assert 1728^2*JJ_1_W12_7 eq formulae_JJ_1;