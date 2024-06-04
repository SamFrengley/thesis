/*
    As promised in the paper we show how we might derive Theorem 3.3 in a highly 
    computational way using the equations of Fisher. 
*/
load "../Z2_3_4/polys.m";

// -------------------------------------
// We first give the direction (Fisher's Equations) => (Thm 3.3)
// -------------------------------------
FF<a,b,u,v> := FunctionField(Rationals(), 4); 

// Here a,b will be so that E1 has c_4(E1) = a and c_6(E1) = b. The [u,v] will be coordinates
// on Fisher's model for X_E^r(N). Our approach is to give alpha and beta which depend only
// on a,b,u,v which solve the polys in Thm 3.3.

//-----------------------
// (N, r) = (2, 1)  -- c.f., Fisher's proof
_, c4, c6 := HessePolynomials(2, 1, [a,b]);
c4 := Evaluate(c4, [u,v]); c6 := Evaluate(c6, [u,v]);

J1 := a^3/(a^3-b^2);
J2 := c4^3/(c4^3-c6^2);

Pa, Pb := CongPolys(2, 1);

alpha := b*(a^3*v^3 - 3*a^2*u^2*v - 3*a*b*u*v^2 - 2*b^2*v^3 - b*u^3)/((3*a*u*v^2 + 2*b*v^3 - u^3)*(a^3 - b^2));
beta := -1*(2)*(a)*(a*u + b*v)*(a^2*v^2 + a*u^2 + 2*b*u*v)/((3*a*u*v^2 +2*b*v^3 - u^3)*(a^3 - b^2));

assert Evaluate(Pa, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta]) eq 0;
assert Evaluate(Pb, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta]) eq 0;

//-----------------------
// (N, r) = (3, 1)
_, c4, c6 := HessePolynomials(3, 1, [a,b]);
c4 := Evaluate(c4, [u,v]); c6 := Evaluate(c6, [u,v]);

J1 := a^3/(a^3-b^2);
J2 := c4^3/(c4^3-c6^2);

Pa, Pb := CongPolys(3, 1);

alpha := 1*(2)*(a)*(3*a^3*v^4 + 2*a*b*u*v^3 - a*u^4 - 2*b^2*v^4 - 2*b*u^3*v)*(3*a^3*v^4 -
        6*a^2*u^2*v^2 - 4*a*b*u*v^3 - a*u^4 - 4*b^2*v^4 - 4*b*u^3*v)/((3*a^2*v^4 +
        6*a*u^2*v^2 + 8*b*u*v^3 - u^4)^(2)*(a^3 - b^2));
beta := 1*(3)*(b)*(a*v^2 - u^2)*(18*a^4*u*v^5 + 9*a^3*b*v^6 + 15*a^2*b*u^2*v^4 +
        6*a^2*u^5*v - 12*a*b^2*u*v^5 + 15*a*b*u^4*v^2 - 8*b^3*v^6 + 20*b^2*u^3*v^3 +
        b*u^6)/((3*a^2*v^4 + 6*a*u^2*v^2 + 8*b*u*v^3 - u^4)^(2)*(a^3 - b^2));
delta := -6*(2*beta^3 - (5*alpha + 2)*beta^2 - 10*(J1-1)*(J2-1)*beta + 3*(J1-1)*(J2-1)*
        (13*alpha - 2 + 6*(J1 + J2)))/(beta^3 - 3*(J1-1)*(J2-1)*beta - 2*((J1-1)*(J2-1))^2);
tau := 36*a^3*v^3 - 36*b^2*v^3;

assert Evaluate(Pa, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta]) eq 0;
assert Evaluate(Pb, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta]) eq 0;
assert tau^2 eq 3*c6*b*delta;

//-----------------------
// (N, r) = (3, 2)
_, c4, c6 := HessePolynomials(3, 2, [a,b]);
c4 := Evaluate(c4, [u,v]); c6 := Evaluate(c6, [u,v]);

J1 := a^3/(a^3-b^2);
J2 := c4^3/(c4^3-c6^2);

Pa, Pb := CongPolys(3, 2);

alpha := 1*(a)*(a^5*v^4 - 6*a^4*u^2*v^2 - 8*a^3*b*u*v^3 - 3*a^3*u^4 - 5*a^2*b^2*v^4 -
        16*a^2*b*u^3*v - 18*a*b^2*u^2*v^2 - 8*b^3*u*v^3 - b^2*u^4)/((a^3 -
        b^2)^(1)*(a^3*v^4 + 6*a^2*u^2*v^2 + 12*a*b*u*v^3 - 3*a*u^4 + 4*b^2*v^4 -
        4*b*u^3*v));
beta := 1*(2)*(3)*(b)*(a*u + b*v)^(2)*(6*a^7*u*v^5 + 5*a^6*b*v^6 + 15*a^5*b*u^2*v^4 +
        18*a^5*u^5*v + 6*a^4*b^2*u*v^5 + 75*a^4*b*u^4*v^2 - 5*a^3*b^3*v^6 +
        140*a^3*b^2*u^3*v^3 + 9*a^3*b*u^6 + 105*a^2*b^3*u^2*v^4 + 30*a^2*b^2*u^5*v +
        36*a*b^4*u*v^5 + 45*a*b^3*u^4*v^2 + 8*b^5*v^6 + 20*b^4*u^3*v^3 - b^3*u^6)/((a^3
        - b^2)^(2)*(a^3*v^4 + 6*a^2*u^2*v^2 + 12*a*b*u*v^3 - 3*a*u^4 + 4*b^2*v^4 -
        4*b*u^3*v)^(2));
tau := -1*(2)^(2)*(3)*(b)*(a*u + b*v)*(6*a^7*u*v^5 + 5*a^6*b*v^6 + 15*a^5*b*u^2*v^4 +
        18*a^5*u^5*v + 6*a^4*b^2*u*v^5 + 75*a^4*b*u^4*v^2 - 5*a^3*b^3*v^6 +
        140*a^3*b^2*u^3*v^3 + 9*a^3*b*u^6 + 105*a^2*b^3*u^2*v^4 + 30*a^2*b^2*u^5*v +
        36*a*b^4*u*v^5 + 45*a*b^3*u^4*v^2 + 8*b^5*v^6 + 20*b^4*u^3*v^3 - b^3*u^6)/((a^3
        - b^2)^(1)*(a^3*v^4 + 6*a^2*u^2*v^2 + 12*a*b*u*v^3 - 3*a*u^4 + 4*b^2*v^4 -
        4*b*u^3*v));

assert Evaluate(Pa, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta]) eq 0;
assert Evaluate(Pb, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta]) eq 0;
assert tau^2 eq 3*b*(2^3*c6)*beta;

//-----------------------
// (N, r) = (4, 1) 
// Note that after doing (4, 1) we get (4, 3) by quadratic twisting by Delta(E)
_, c4, c6 := HessePolynomials(4, 1, [a,b]);
c4 := Evaluate(c4, [u,v]); c6 := Evaluate(c6, [u,v]);

J1 := a^3/(a^3-b^2);
J2 := c4^3/(c4^3-c6^2);

Pa, Pb, Pg := CongPolys(4, 1);

alpha := -1*(b)*(2916*a^7*u*v^11 + 1215*a^6*b*v^12 + 1188*a^6*u^3*v^9 +
        4158*a^5*b*u^2*v^10 + 2376*a^5*u^5*v^7 - 4464*a^4*b^2*u*v^11 +
        5841*a^4*b*u^4*v^8 - 792*a^4*u^7*v^5 - 2240*a^3*b^3*v^12 + 4224*a^3*b^2*u^3*v^9
        - 924*a^3*b*u^6*v^6 - 44*a^3*u^9*v^3 - 4224*a^2*b^3*u^2*v^10 -
        3168*a^2*b^2*u^5*v^7 - 495*a^2*b*u^8*v^4 - 12*a^2*u^11*v + 1536*a*b^4*u*v^11 -
        6336*a*b^3*u^4*v^8 - 66*a*b*u^10*v^2 + 1024*b^5*v^12 - 5632*b^4*u^3*v^9 -
        176*b^2*u^9*v^3 - b*u^12)/((a^3 - b^2)^(1)*(27*a^3*v^6 - 45*a^2*u^2*v^4 -
        24*a*b*u*v^5 - 15*a*u^4*v^2 - 32*b^2*v^6 - 40*b*u^3*v^3 + u^6)^(2));
beta := 1*(2)*(a)*(9*a^3*v^4 + 6*a^2*u^2*v^2 + 4*a*b*u*v^3 + a*u^4 - 8*b^2*v^4 +
        4*b*u^3*v)*(81*a^5*v^8 + 252*a^4*u^2*v^6 + 264*a^3*b*u*v^7 - 42*a^3*u^4*v^4 -
        80*a^2*b^2*v^8 + 56*a^2*b*u^3*v^5 + 28*a^2*u^6*v^2 - 224*a*b^2*u^2*v^6 +
        56*a*b*u^5*v^3 + a*u^8 - 256*b^3*u*v^7 + 112*b^2*u^4*v^4 + 8*b*u^7*v)/((a^3 -
        b^2)^(1)*(27*a^3*v^6 - 45*a^2*u^2*v^4 - 24*a*b*u*v^5 - 15*a*u^4*v^2 - 32*b^2*v^6
        - 40*b*u^3*v^3 + u^6)^(2));
gamma := 1*(2)^(2)*(3)*(v)^(2)*(a)^(2)*(81*a^5*v^8 + 252*a^4*u^2*v^6 + 264*a^3*b*u*v^7 -
        42*a^3*u^4*v^4 - 80*a^2*b^2*v^8 + 56*a^2*b*u^3*v^5 + 28*a^2*u^6*v^2 -
        224*a*b^2*u^2*v^6 + 56*a*b*u^5*v^3 + a*u^8 - 256*b^3*u*v^7 + 112*b^2*u^4*v^4 +
        8*b*u^7*v)^(2)/((a^3 - b^2)^(1)*(27*a^3*v^6 - 45*a^2*u^2*v^4 - 24*a*b*u*v^5 -
        15*a*u^4*v^2 - 32*b^2*v^6 - 40*b*u^3*v^3 + u^6)^(3));
delta := 3*(J1*J2 - (alpha + 1)^2);
tau := -1*(2)^(2)*(3)*(v)*(b)*(3*a*u*v^2 + 2*b*v^3 - u^3)*(243*a^6*v^8 +
        180*a^5*u^2*v^6 + 24*a^4*b*u*v^7 + 258*a^4*u^4*v^4 - 496*a^3*b^2*v^8 +
        424*a^3*b*u^3*v^5 + 20*a^3*u^6*v^2 - 96*a^2*b^2*u^2*v^6 + 168*a^2*b*u^5*v^3 +
        3*a^2*u^8 - 48*a*b^2*u^4*v^4 + 24*a*b*u^7*v + 256*b^4*v^8 - 256*b^3*u^3*v^5 +
        64*b^2*u^6*v^2)*(2916*a^7*u*v^11 + 1215*a^6*b*v^12 + 1188*a^6*u^3*v^9 +
        4158*a^5*b*u^2*v^10 + 2376*a^5*u^5*v^7 - 4464*a^4*b^2*u*v^11 +
        5841*a^4*b*u^4*v^8 - 792*a^4*u^7*v^5 - 2240*a^3*b^3*v^12 + 4224*a^3*b^2*u^3*v^9
        - 924*a^3*b*u^6*v^6 - 44*a^3*u^9*v^3 - 4224*a^2*b^3*u^2*v^10 -
        3168*a^2*b^2*u^5*v^7 - 495*a^2*b*u^8*v^4 - 12*a^2*u^11*v + 1536*a*b^4*u*v^11 -
        6336*a*b^3*u^4*v^8 - 66*a*b*u^10*v^2 + 1024*b^5*v^12 - 5632*b^4*u^3*v^9 -
        176*b^2*u^9*v^3 - b*u^12)/((a^3 - b^2)^(1)*(27*a^3*v^6 - 45*a^2*u^2*v^4 -
        24*a*b*u*v^5 - 15*a*u^4*v^2 - 32*b^2*v^6 - 40*b*u^3*v^3 + u^6)^(3));

assert Evaluate(Pa, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta, gamma]) eq 0;
assert Evaluate(Pb, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta, gamma]) eq 0;
assert Evaluate(Pg, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta, gamma]) eq 0;
assert tau^2 eq 3*b*c6*alpha*delta;

// -------------------------------------
// We now give the direction (Thm 3.3) => (Fisher's equations)
// -------------------------------------
KK<J1,J2> := FunctionField(Rationals(), 2);
cc := [1728*J1/(1728*J1-1728), 1728*J1/(1728*J1-1728)];

// We will start with the quotient ring K(J1,J2)[alpha,beta,gamma]/(Pa,Pb,Pg). Then cook up coordinates
// [u,v] on X_E^r(N) = P^1.

//-----------------------
// (N, r) = (2, 1)  -- c.f., Fisher's proof
Pa, Pb := CongPolys(2, 1);

PP<alpha,beta> := PolynomialRing(KK,2);
LL<alpha,beta> := quo<PP| Evaluate(Pa, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta]),
                          Evaluate(Pb, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta])>;
u := (beta + J1)*(1 - alpha - J1);
v := (J1 - J2)*(J1-1);
_,c4,c6 := HessePolynomials(2, 1, cc : Variables := [u,v]);
E2 := EllipticCurve([-27*c4,-54*c6]);
assert jInvariant(E2) eq 1728*J2;

//-----------------------
// (N, r) = (4, 1) 
Pa, Pb, Pg := CongPolys(4, 1);

// In this case I couldn't set up the equations to work quickly enough, so using 
// the parametrisation in Thm 3.3. 
KJ<J1> := FunctionField(Rationals());
cc := [1728*J1/(1728*J1-1728), 1728*J1/(1728*J1-1728)];

A2_K<t_2,g> := AffineSpace(KJ, 2);              // t_2 is u/v in the (2, 1) case above, g is gamma
A1_K<t_4> := AffineSpace(KJ, 1);                // t_4 is playing the role of u/v

// The (2, 1) parametrisation
J2 := 1*(J1)*(J1*t_2^2 + 2*J1*t_2 + J1 - t_2^2 - 2*t_2)^(3)/((J1 - 1)^(1)*(J1*t_2^3 - 3*J1*t_2 - 
2*J1 - t_2^3)^(2));
a := 1*(J1^2*t_2^3 + 3*J1^2*t_2^2 + 3*J1^2*t_2 + J1^2 - 2*J1*t_2^3 - 3*J1*t_2^2 - 3*J1*t_2 - 2*J1
+ t_2^3)/((J1*t_2^3 - 3*J1*t_2 - 2*J1 - t_2^3)^(1));
b := 1*(2)*(t_2 + 1)*(J1)*(J1*t_2^2 + 2*J1*t_2 + J1 - t_2^2 - 2*t_2)/((J1*t_2^3 - 3*J1*t_2 - 2*J1 -
t_2^3)^(1));

Pa, Pb, Pg := CongPolys(4, 1);
assert Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), a, b, g]) eq 0;
assert Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), a, b, g]) eq 0;
Pg := Evaluate(Pg, [J1*J2, (J1-1)*(J2-1), a, b, g]);

calZ41 := Curve(A2_K, Numerator(Pg));

// CLAIM: This parametrisation actually gives us Fisher's Hesse family
f1 := map<calZ41 -> A2_K | [t_2, 
            ((J1 - 1)*(J1*t_2^3 - 3*J1*t_2 - 2*J1 - t_2^3)^2/((J1*t_2^2 + 2*J1*t_2 + J1 - t_2^2 - 2*t_2)^2*J1))*g]>;
X := Image(f1);
f1 := iso<calZ41 -> X | DefiningPolynomials(f1), 
                        [t_2, X.2/((J1 - 1)*(J1*t_2^3 - 3*J1*t_2 - 2*J1 - t_2^3)^2/((J1*t_2^2 + 2*J1*t_2 + J1 - t_2^2 - 2*t_2)^2*J1))]>;

f2 := (16*J1^5*t_2^7 + 24*J1^5*t_2^6 - 66*J1^5*t_2^5 - 148*J1^5*t_2^4 - 12*J1^5*t_2^3 + 
        168*J1^5*t_2^2 + 142*J1^5*t_2 + 36*J1^5 - 64*J1^4*t_2^7 - 72*J1^4*t_2^6 + 12*J1^4*t_2^5*(X.2) 
        + 198*J1^4*t_2^5 + 24*J1^4*t_2^4*(X.2) + 340*J1^4*t_2^4 - 24*J1^4*t_2^3*(X.2) + 24*J1^4*t_2^3 - 
        96*J1^4*t_2^2*(X.2) - 240*J1^4*t_2^2 - 84*J1^4*t_2*(X.2) - 166*J1^4*t_2 - 24*J1^4*(X.2) - 36*J1^4 + 
        96*J1^3*t_2^7 + 72*J1^3*t_2^6 - 36*J1^3*t_2^5*(X.2) - 198*J1^3*t_2^5 - 58*J1^3*t_2^4*(X.2) - 
        236*J1^3*t_2^4 + 48*J1^3*t_2^3*(X.2) - 12*J1^3*t_2^3 + 150*J1^3*t_2^2*(X.2) + 72*J1^3*t_2^2 + 
        104*J1^3*t_2*(X.2) + 24*J1^3*t_2 + 24*J1^3*(X.2) - 64*J1^2*t_2^7 - 24*J1^2*t_2^6 + 36*J1^2*t_2^5*(X.2)
        + 66*J1^2*t_2^5 + 44*J1^2*t_2^4*(X.2) + 44*J1^2*t_2^4 - 24*J1^2*t_2^3*(X.2) - J1^2*t_2^2*(X.2)^2 - 
        54*J1^2*t_2^2*(X.2) - 2*J1^2*t_2*(X.2)^2 - 20*J1^2*t_2*(X.2) - J1^2*(X.2)^2 + 16*J1*t_2^7 - 12*J1*t_2^5*(X.2) 
        - 10*J1*t_2^4*(X.2) + J1*t_2^2*(X.2)^2 - J1*t_2*(X.2)^3 + 2*J1*t_2*(X.2)^2 - J1*(X.2)^3 + 
        t_2*(X.2)^3)/((2)^(1)*(J1)*(J1 - 1)*(J1*t_2^3 - 3*J1*t_2 - 2*J1 - t_2^3)*(2*J1^2*t_2^3 + 
        3*J1^2*t_2^2 - J1^2 - 4*J1*t_2^3 - 3*J1*t_2^2 + 2*t_2^3));
f2_inv := [(J1^2*t_4^4 + 6*J1^2*t_4^2 + 16*J1^2*t_4 + 9*J1^2 - 2*J1*t_4^4 - 6*J1*t_4^2 - 16*J1*t_4 + 
        t_4^4)/((2)^(2)*(J1 - 1)*(J1*t_4^3 - 3*J1*t_4 - 2*J1 - t_4^3)),
        (3)*(J1)*(J1^3*t_4^6 - 15*J1^3*t_4^4 - 40*J1^3*t_4^3 - 45*J1^3*t_4^2 - 24*J1^3*t_4 - 
        5*J1^3 - 3*J1^2*t_4^6 + 30*J1^2*t_4^4 + 80*J1^2*t_4^3 + 45*J1^2*t_4^2 + 24*J1^2*t_4 + 
        32*J1^2 + 3*J1*t_4^6 - 15*J1*t_4^4 - 40*J1*t_4^3 - t_4^6)/((2)^(2)*(J1 - 1)*(J1*t_4^3 - 
        3*J1*t_4 - 2*J1 - t_4^3)^(2))];
f2 := iso<X -> A1_K | [f2], f2_inv>; 

// To clarify: given a (J1, J2, alpha, beta, gamma) which satisfy the equations we can write
// t_2 := (beta + J1)*(1 - alpha - J1)/((J1 - J2)*(J1-1)); (from the (2, 1) calculation)
// then putting t_4 to be the image of (t_2, gamma) under f1*f2 gives the "t_4" which would
// be put into HessePolynomials(4, 1, cc : variables:=[t,1]). 
// To prove this note that we have:

_<t_4> := FunctionField(Rationals());
t_2 := DefiningPolynomials(Inverse(f2)*Inverse(f1))[1];
t_2 := Evaluate(t_2, [t_4]);
J2 := 1*(J1)*(J1*t_2^2 + 2*J1*t_2 + J1 - t_2^2 - 2*t_2)^(3)/((J1 - 1)^(1)*(J1*t_2^3 - 3*J1*t_2 - 
2*J1 - t_2^3)^(2));

_,c4,c6 := HessePolynomials(4, 1, cc : Variables := [t_4,1]);
E2 := EllipticCurve([-27*c4,-54*c6]);
assert jInvariant(E2) eq 1728*J2;


/*
// You can print out the formulae for t_4 in terms of t_2 and gamma with the following commands:

pp := DefiningPolynomials(f1*f2);
_<t_2, gamma> := FunctionField(KJ, 4);
t_4 := Evaluate(pp[1], [t_2, gamma]);

*/

//-----------------------
// (N, r) = (3, 1)
KK<J1,J2> := FunctionField(Rationals(), 2);
cc := [1728*J1/(1728*J1-1728), 1728*J1/(1728*J1-1728)];

Pa, Pb := CongPolys(3, 1);

PP<alpha,beta> := PolynomialRing(KK,2);
LL<alpha,beta> := quo<PP| Evaluate(Pa, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta]),
                          Evaluate(Pb, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta])>;
u := -1*(J1)*(8*J1^6*J2^2 - 16*J1^6*J2 + 8*J1^6 + 32*J1^5*J2^3 + 24*J1^5*J2^2*alpha - 
        75*J1^5*J2^2 - 48*J1^5*J2*alpha + 54*J1^5*J2 + 24*J1^5*alpha - 11*J1^5 + 104*J1^4*J2^4 +
        120*J1^4*J2^3*alpha - 258*J1^4*J2^3 + 12*J1^4*J2^2*alpha^2 - 276*J1^4*J2^2*alpha - 
        8*J1^4*J2^2*beta^2 + 198*J1^4*J2^2 - 24*J1^4*J2*alpha^2 - 4*J1^4*J2*alpha*beta^2 + 
        192*J1^4*J2*alpha + 18*J1^4*J2*beta^2 - 38*J1^4*J2 + 12*J1^4*alpha^2 + 4*J1^4*alpha*beta^2 - 
        36*J1^4*alpha - 10*J1^4*beta^2 - 6*J1^4 - 227*J1^3*J2^4 - 120*J1^3*J2^3*alpha^2 - 
        252*J1^3*J2^3*alpha + 8*J1^3*J2^3*beta^2 + 466*J1^3*J2^3 - 72*J1^3*J2^2*alpha^3 + 
        204*J1^3*J2^2*alpha^2 + 4*J1^3*J2^2*alpha*beta^2 + 504*J1^3*J2^2*alpha - 12*J1^3*J2^2*beta^2 - 
        238*J1^3*J2^2 + 144*J1^3*J2*alpha^3 + 4*J1^3*J2*alpha^2*beta^2 - 48*J1^3*J2*alpha^2 + 
        4*J1^3*J2*alpha*beta^2 - 252*J1^3*J2*alpha + 2*J1^3*J2*beta^2 - 14*J1^3*J2 - 72*J1^3*alpha^3 - 
        4*J1^3*alpha^2*beta^2 - 36*J1^3*alpha^2 - 8*J1^3*alpha*beta^2 + J1^3*beta^4 + 2*J1^3*beta^2 + 13*J1^3 + 
        138*J1^2*J2^4 + 252*J1^2*J2^3*alpha^2 + 144*J1^2*J2^3*alpha - 6*J1^2*J2^3*beta^2 - 
        278*J1^2*J2^3 + 36*J1^2*J2^2*alpha^4 + 144*J1^2*J2^2*alpha^3 - 4*J1^2*J2^2*alpha^2*beta^2 - 
        468*J1^2*J2^2*alpha^2 - 8*J1^2*J2^2*alpha*beta^2 - 276*J1^2*J2^2*alpha + 2*J1^2*J2^2*beta^2 + 
        138*J1^2*J2^2 - 72*J1^2*J2*alpha^4 - 288*J1^2*J2*alpha^3 - 4*J1^2*J2*alpha^2*beta^2 + 
        180*J1^2*J2*alpha^2 + 4*J1^2*J2*alpha*beta^2 + 120*J1^2*J2*alpha - 2*J1^2*J2*beta^4 - 
        4*J1^2*J2*beta^2 + 6*J1^2*J2 + 36*J1^2*alpha^4 + 144*J1^2*alpha^3 + 8*J1^2*alpha^2*beta^2 + 
        36*J1^2*alpha^2 + 4*J1^2*alpha*beta^2 + 12*J1^2*alpha - 4*J1^2*beta^4 + 8*J1^2*beta^2 - 4*J1^2 - 
        11*J1*J2^4 - 144*J1*J2^3*alpha^2 - 12*J1*J2^3*alpha + 6*J1*J2^3*beta^2 + 30*J1*J2^3 - 
        72*J1*J2^2*alpha^4 - 72*J1*J2^2*alpha^3 + 8*J1*J2^2*alpha^2*beta^2 + 276*J1*J2^2*alpha^2 + 
        4*J1*J2^2*alpha*beta^2 + 24*J1*J2^2*alpha + J1*J2^2*beta^4 + 10*J1*J2^2*beta^2 - 27*J1*J2^2 + 
        144*J1*J2*alpha^4 + 144*J1*J2*alpha^3 - 4*J1*J2*alpha^2*beta^2 - 120*J1*J2*alpha^2 - 4*J1*J2*alpha*beta^2 
        - 12*J1*J2*alpha + 8*J1*J2*beta^4 - 16*J1*J2*beta^2 + 8*J1*J2 - 72*J1*alpha^4 - 72*J1*alpha^3 - 
        4*J1*alpha^2*beta^2 - 12*J1*alpha^2 - 4*J2^4 + 12*J2^3*alpha^2 - 8*J2^3*beta^2 + 8*J2^3 + 
        36*J2^2*alpha^4 - 4*J2^2*alpha^2*beta^2 - 24*J2^2*alpha^2 - 4*J2^2*beta^4 + 8*J2^2*beta^2 - 4*J2^2 - 
        72*J2*alpha^4 + 4*J2*alpha^2*beta^2 + 12*J2*alpha^2 + 36*alpha^4);
v := ((2)^(1)*(J1 - 1)*(J1 - J2)*(4*J1^5*J2^2 - 8*J1^5*J2 + 4*J1^5 + 24*J1^4*J2^3 + 14*J1^4*J2^2*alpha - 
        51*J1^4*J2^2 - 28*J1^4*J2*alpha + 30*J1^4*J2 + 14*J1^4*alpha - 3*J1^4 - 45*J1^3*J2^3 - 
        14*J1^3*J2^2*alpha^2 - 24*J1^3*J2^2*alpha - 4*J1^3*J2^2*beta^2 + 84*J1^3*J2^2 + 
        28*J1^3*J2*alpha^2 - 2*J1^3*J2*alpha*beta^2 + 48*J1^3*J2*alpha - 2*J1^3*J2*beta^2 - 33*J1^3*J2 - 
        14*J1^3*alpha^2 + 2*J1^3*alpha*beta^2 - 24*J1^3*alpha + 6*J1^3*beta^2 - 6*J1^3 + 18*J1^2*J2^3 + 
        24*J1^2*J2^2*alpha^2 + 6*J1^2*J2^2*alpha + 2*J1^2*J2^2*beta^2 - 31*J1^2*J2^2 + 
        2*J1^2*J2*alpha^2*beta^2 - 48*J1^2*J2*alpha^2 - 2*J1^2*J2*alpha*beta^2 - 12*J1^2*J2*alpha + 
        4*J1^2*J2*beta^2 + 8*J1^2*J2 - 2*J1^2*alpha^2*beta^2 + 24*J1^2*alpha^2 + 2*J1^2*alpha*beta^2 + 
        6*J1^2*alpha + J1^2*beta^4 - 6*J1^2*beta^2 + 5*J1^2 + 3*J1*J2^3 - 6*J1*J2^2*alpha^2 + 
        4*J1*J2^2*alpha + 2*J1*J2^2*beta^2 - 6*J1*J2^2 + 2*J1*J2*alpha^2*beta^2 + 12*J1*J2*alpha^2 + 
        4*J1*J2*alpha*beta^2 - 8*J1*J2*alpha - J1*J2*beta^4 - 2*J1*J2*beta^2 + 3*J1*J2 - 2*J1*alpha^2*beta^2 - 
        6*J1*alpha^2 - 4*J1*alpha*beta^2 + 4*J1*alpha - 4*J2^2*alpha^2 - 4*J2*alpha^2*beta^2 + 8*J2*alpha^2 + 
        4*alpha^2*beta^2 - 4*alpha^2));
_,c4,c6 := HessePolynomials(3, 1, cc : Variables := [u,v]);
E2 := EllipticCurve([-27*c4,-54*c6]);
assert jInvariant(E2) eq 1728*J2;

//-----------------------
// (N, r) = (3, 2)
Pa, Pb := CongPolys(3, 2);

PP<alpha,beta> := PolynomialRing(KK,2);
LL<alpha,beta> := quo<PP| Evaluate(Pa, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta]),
                          Evaluate(Pb, [J1*J2, (J1 - 1)*(J2 - 1), alpha, beta])>;
u := -1*(16*J1^4*J2^2 - 32*J1^4*J2 + 16*J1^4 - 9*J1^3*J2^2*alpha^2 + 
        2*J1^3*J2^2*alpha - 57*J1^3*J2^2 + 18*J1^3*J2*alpha^2 - 4*J1^3*J2*alpha + 
        114*J1^3*J2 - 9*J1^3*alpha^2 + 2*J1^3*alpha - 57*J1^3 + 18*J1^2*J2^2*alpha^2 + 
        12*J1^2*J2^2*alpha + 66*J1^2*J2^2 - 36*J1^2*J2*alpha^2 + 10*J1^2*J2*alpha*beta^2
        - 24*J1^2*J2*alpha - 6*J1^2*J2*beta^2 - 132*J1^2*J2 + 18*J1^2*alpha^2 - 
        10*J1^2*alpha*beta^2 + 12*J1^2*alpha + 6*J1^2*beta^2 + 66*J1^2 - 
        9*J1*J2^2*alpha^2 - 30*J1*J2^2*alpha - 25*J1*J2^2 + 18*J1*J2*alpha^2 - 
        26*J1*J2*alpha*beta^2 + 60*J1*J2*alpha + 6*J1*J2*beta^2 + 50*J1*J2 - 
        9*J1*alpha^2 + 26*J1*alpha*beta^2 - 30*J1*alpha - J1*beta^4 - 6*J1*beta^2 - 
        25*J1 + 16*J2^2*alpha + 16*J2*alpha*beta^2 - 32*J2*alpha - 16*alpha*beta^2 + 
        16*alpha);
v := ((16*J1^4*J2^2 - 32*J1^4*J2 + 16*J1^4 - 9*J1^3*J2^2*alpha^2 + 
        6*J1^3*J2^2*alpha - 61*J1^3*J2^2 + 18*J1^3*J2*alpha^2 - 12*J1^3*J2*alpha + 
        122*J1^3*J2 - 9*J1^3*alpha^2 + 6*J1^3*alpha - 61*J1^3 + 18*J1^2*J2^2*alpha^2 + 
        78*J1^2*J2^2 - 36*J1^2*J2*alpha^2 + 10*J1^2*J2*alpha*beta^2 - 10*J1^2*J2*beta^2 
        - 156*J1^2*J2 + 18*J1^2*alpha^2 - 10*J1^2*alpha*beta^2 + 10*J1^2*beta^2 + 
        78*J1^2 - 9*J1*J2^2*alpha^2 - 18*J1*J2^2*alpha - 37*J1*J2^2 + 18*J1*J2*alpha^2 -
        38*J1*J2*alpha*beta^2 + 36*J1*J2*alpha + 2*J1*J2*beta^2 + 74*J1*J2 - 
        9*J1*alpha^2 + 38*J1*alpha*beta^2 - 18*J1*alpha - J1*beta^4 - 2*J1*beta^2 - 
        37*J1 + 12*J2^2*alpha + 4*J2^2 + 28*J2*alpha*beta^2 - 24*J2*alpha + 8*J2*beta^2 
        - 8*J2 - 28*alpha*beta^2 + 12*alpha + 4*beta^4 - 8*beta^2 + 4)^(1));
_,c4,c6 := HessePolynomials(3, 2, cc : Variables := [u,v]);
E2 := EllipticCurve([-27*2^2*c4,-54*2^3*c6]);
assert jInvariant(E2) eq 1728*J2;
