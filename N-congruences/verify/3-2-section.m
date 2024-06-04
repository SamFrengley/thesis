/*
    Computational verifications for Section 5.3.
*/
load "../Z12-r/Z12equations.m";
load "../Z2_3_4/polys.m";   

QQ := Rationals();

Pa_3, Pb_3 := CongPolys(3, 2);
Pa_4, Pb_4, Pg_4 := CongPolys(4, 1);

//---------------------------
// W(4, r)
W4r<a,b> := AffineSpace(QQ, 2);
_<x> := PolynomialRing(FieldOfFractions(CoordinateRing(W4r)));

JJ := -1*(2)^(4)*(2*a*b - 2*a - b^2 - 2*b)^(3)/((3)^(3)*(b - 1)^(2)*(a - b - 1)^(4));
JJ_1 := 1*(a^2*b - a^2 + 2*a*b^2 - 2*a + b^3 + b^2 + b + 1)^(2)/((b - 1)^(2)*(a - b - 1)^(4));

A5<alpha_4_var, beta_4_var, gamma_4_var, JJ_var, JJ_1_var> := AffineSpace(QQ, 5);

calW4r := Surface(A5, [ Evaluate(Pa_4, [JJ_var, JJ_1_var, alpha_4_var, beta_4_var, gamma_4_var]), 
                        Evaluate(Pb_4, [JJ_var, JJ_1_var, alpha_4_var, beta_4_var, gamma_4_var]), 
                        Evaluate(Pg_4, [JJ_var, JJ_1_var, alpha_4_var, beta_4_var, gamma_4_var])]);


alpha_4 := Roots(Evaluate(Pa_4, [JJ, JJ_1, x, 0, 0]))[1][1];
beta_4 := Roots(Evaluate(Pb_4, [JJ, JJ_1, alpha_4, x, 0]))[1][1];
gamma_4 := Roots(Evaluate(Pg_4, [JJ, JJ_1, alpha_4, beta_4, x]))[1][1];

phi := iso<W4r -> calW4r | [alpha_4, beta_4, gamma_4, JJ, JJ_1], [(2*beta_4_var*gamma_4_var*JJ_var + 4*JJ_var^2)/(beta_4_var*gamma_4_var*JJ_var + gamma_4_var^3), (beta_4_var*JJ_var - gamma_4_var^2)/(beta_4_var*JJ_var + gamma_4_var^2)] : Check:=true, CheckInverse:=true>;

//--------------------------
// W(4, r)^{cbrt(DD')}
A3<a,b,alpha_3> := AffineSpace(QQ, 3); 
S1 := Surface(A3, Numerator(Numerator(Evaluate(Pa_3, [Evaluate(JJ, [a,b]), Evaluate(JJ_1, [a,b]), alpha_3, 0]))));

A3<a,b,alpha_3_prime> := AffineSpace(QQ, 3); 
S2 := Surface(A3, alpha_3_prime^3 - 2*(b - 1)*(-a + b + 1)^2);

alpha_3 := 2*(-2*a*b + b^2 + 2*a + 2*b)/(3*(b-1)*(-a + b + 1)^2)*alpha_3_prime;
phi := map<S2 -> S1 | [a,b,alpha_3]>;
assert IsInvertible(phi);

W4r_cbrtDD<c,d> := AffineSpace(QQ, 2);

a := (-1/2*c*d^3 + 4*d^3 + 2)/(c + 2*d^3 + 1);
b := (-c + 2*d^3 + 1)/(c + 2*d^3 + 1);
alpha_3_prime := -c*d^2/(c + 2*d^3 + 1);

phi := map<W4r_cbrtDD -> S2 | [a, b, alpha_3_prime]>;
assert IsInvertible(phi);

JJ := Evaluate(JJ, [a, b]);
JJ_1 := Evaluate(JJ_1, [a, b]);
alpha_4 := Evaluate(alpha_4, [a, b]);
alpha_3 := Evaluate(alpha_3, [a, b, alpha_3_prime]);

//--------------------------
// W^+(12, 5)
A4<c,d,beta_3,theta> := AffineSpace(QQ, 4); 

f1 := beta_3^2 - theta;
f2 := theta^2 -6*Evaluate((alpha_3 + 1)*(JJ_1), [c,d])*theta -8*Evaluate((JJ_1)^2, [c,d])*beta_3 - 3*Evaluate((alpha_3 - 1)^2*(JJ_1)^2, [c,d]);

A4<c,d,beta_3_prime,theta_prime> := AffineSpace(QQ, 4); 

g1 := c^2*d^4 - 8*c*d^5 - 4*c*d^2 - 8*d^3*theta_prime - 8*d^3*beta_3_prime - theta_prime^2 - 2*theta_prime*beta_3_prime - 4*theta_prime + 2*beta_3_prime^2 - 4*beta_3_prime;
g2 := -3*c^2*d^4 + 8*c^2*d^3 + 4*c^2 - 48*c*d^3 - 6*c*d^2*beta_3_prime - 24*c + 3*theta_prime^2 - 3*beta_3_prime^2;

beta_3 := (8*c^2*d^3 + 4*c^2 - 12*c*d^5 - 48*c*d^3 - 3*c*d^2*beta_3_prime - 6*c*d^2 - 24*c + 48*d^6
        - 24*d^3*beta_3_prime + 48*d^3 - 12*beta_3_prime + 12)*(c^3*d^6 - 32*c^2*d^6 - 32*c^2*d^3 - 8*c^2 + 
        192*c*d^6 + 192*c*d^3 + 48*c - 64*d^9 - 96*d^6 - 48*d^3 - 
        8)/((3)^(1)*(d)^(8)*(4*d^3 + beta_3_prime + 2)*(c)^(4));
theta := (12*c^2*d^7 - 32*c^2*d^6 + 3*c^2*d^4*beta_3_prime + 6*c^2*d^4 - 8*c^2*d^3*beta_3_prime + 
        16*c^2*d^3*theta_prime - 32*c^2*d^3 - 4*c^2*beta_3_prime + 8*c^2*theta_prime - 8*c^2 + 192*c*d^8 + 192*c*d^6 + 
        48*c*d^5*beta_3_prime + 192*c*d^5 + 48*c*d^3*beta_3_prime - 96*c*d^3*theta_prime + 192*c*d^3 + 24*c*d^2*beta_3_prime + 
        48*c*d^2 + 24*c*beta_3_prime - 48*c*theta_prime + 48*c + 192*d^9 + 48*d^6*beta_3_prime + 288*d^6*theta_prime + 288*d^6 + 
        48*d^3*beta_3_prime + 288*d^3*theta_prime + 144*d^3 + 12*beta_3_prime + 72*theta_prime + 24)*(c^3*d^6 - 32*c^2*d^6 - 
        32*c^2*d^3 - 8*c^2 + 192*c*d^6 + 192*c*d^3 + 48*c - 64*d^9 - 96*d^6 - 48*d^3 - 
        8)^(2)/((3)^(1)*(d)^(16)*(4*d^3 + beta_3_prime + 2)*(c)^(8));

assert Numerator(Evaluate(f1, [c, d, beta_3, theta])) in Ideal([g1, g2]);
assert Numerator(Evaluate(f2, [c, d, beta_3, theta])) in Ideal([g1, g2]);

S := Surface(A4, [g1, g2]);
i1 := -(3*c*d^3 - 5*c*d^2 + 4*c*d - 2*c + 3*d*beta_3_prime + 3*d*theta_prime + 3*beta_3_prime - 3*theta_prime)/((2*c*d^2 - 4*c*d + 2*c + 3*beta_3_prime));
i2 := 3*(d)*(c*d^2 - 2*c*d + 2*c + beta_3_prime + theta_prime)/((2*c*d^2 - 4*c*d + 2*c + 3*beta_3_prime));

W12_5_plus<s,t> := AffineSpace(QQ, 2);

c := 6*(s + t)*(s^3 + 3*s^2*t + 6*s^2 + 3*s*t^2 + 12*s*t + 12*s + 3*t^3 + 
    6*t^2 + 12*t + 8)/((2*s^4 + 4*s^3*t + 4*s^3 + 5*s^2*t^2 + 18*s^2*t + 12*s^2 + 
    6*s*t^3 + 20*s*t^2 + 24*s*t + 10*s + 3*t^4 + 6*t^3 + 11*t^2 + 8*t - 1));
d := t/(s + t + 2);
beta_3_prime := -1*(2)^(2)*(s - 1)*(s^2 + s*t + s + 2*t + 1)*(s^3 + 3*s^2*t + 6*s^2 + 3*s*t^2 + 
    12*s*t + 12*s + 3*t^3 + 6*t^2 + 12*t + 8)/((s + t + 2)^(2)*(2*s^4 + 4*s^3*t + 
    4*s^3 + 5*s^2*t^2 + 18*s^2*t + 12*s^2 + 6*s*t^3 + 20*s*t^2 + 24*s*t + 10*s + 
    3*t^4 + 6*t^3 + 11*t^2 + 8*t - 1));
theta_prime := 1*(2)*(s^3 + 3*s^2*t + 6*s^2 + 3*s*t^2 + 12*s*t + 12*s + 3*t^3 + 6*t^2 + 12*t + 
    8)*(4*s^3 + 4*s^2*t - 3*s*t^2 - 2*s*t - 6*s - 3*t^3 - 2*t + 2)/((s + t + 
    2)^(2)*(2*s^4 + 4*s^3*t + 4*s^3 + 5*s^2*t^2 + 18*s^2*t + 12*s^2 + 6*s*t^3 + 
    20*s*t^2 + 24*s*t + 10*s + 3*t^4 + 6*t^3 + 11*t^2 + 8*t - 1));

phi := iso<W12_5_plus -> S | [c, d, beta_3_prime, theta_prime], [i1, i2] : Check:=true, CheckInverse:=true>;

JJ := Evaluate(JJ, [c, d]);
JJ_1 := Evaluate(JJ_1, [c, d]);
alpha_4 := Evaluate(alpha_4, [c, d]);
alpha_3 := Evaluate(alpha_3, [c, d]);
beta_3 := Evaluate(beta_3, [c, d, beta_3_prime, theta_prime]);
JpJ := JJ - JJ_1 + 1;
delta_4 := 3*(JJ - (alpha_4 + 1)^2);

//--------------------------
// W(12, 5)
A3<s,t,w> := AffineSpace(QQ, 3); 

f := -s^2 + t^2 + 1;
assert IsSquare(Evaluate(alpha_4*beta_3*delta_4, [s,t])*f);
S := Surface(A3, w^2 - f);

W12_5<u,v> := AffineSpace(QQ, 2);

s := (-u^2 + 2*u + v^2 - 2)/(u^2 - 2*u - v^2);
t := (2*u - 2)/(u^2 - 2*u - v^2);
w := 2*v/(u^2 - 2*u - v^2);

phi := map<W12_5 -> S | [s, t, w]>;
assert IsInvertible(phi);

JJ_W12_5 := Evaluate(JJ, [s, t]);
JJ_1_W12_5 := Evaluate(JJ_1, [s, t]);
JpJ_W12_5 := JJ_W12_5 - JJ_1_W12_5 + 1;

assert IsSquare((JpJ_W12_5^2 - 4*JJ_W12_5)*Z12Equations(5, u, v));

formulae_JJ, formulae_JJ_1 := W12Moduli(5, u, v);
assert 1728^2*JJ_W12_5 eq formulae_JJ;
assert 1728^2*JJ_1_W12_5 eq formulae_JJ_1;

//--------------------------
// W(12, 11)
A3<s,t,w> := AffineSpace(QQ, 3); 

f := -(s + t)*(s + t + 2);
assert IsSquare(Evaluate(alpha_4*beta_3, [s,t])*f);
S := Surface(A3, w^2 - f);

W12_11<u,v> := AffineSpace(QQ, 2);

s := (-u*v^2 + u - 3*v^2 - 1)/(u*v^2 + u + v^2 + 1);
t := (-u + 1)/(u + 1);
w := 2*v/(v^2 + 1);

phi := map<W12_11 -> S | [s, t, w]>;
assert IsInvertible(phi);

JJ_W12_11 := Evaluate(JJ, [s, t]);
JJ_1_W12_11 := Evaluate(JJ_1, [s, t]);
JpJ_W12_11 := JJ_W12_11 - JJ_1_W12_11 + 1;

assert IsSquare((JpJ_W12_11^2 - 4*JJ_W12_11)*Z12Equations(11, u, v));

formulae_JJ, formulae_JJ_1 := W12Moduli(11, u, v);
assert 1728^2*JJ_W12_11 eq formulae_JJ;
assert 1728^2*JJ_1_W12_11 eq formulae_JJ_1;