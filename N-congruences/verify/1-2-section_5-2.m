/*
    Computational verifications for Section 3.2, when (N,r) = (5,2).
*/
QQ := Rationals();
KK<z> := CyclotomicField(5);
PP<x0,x1,y0,y1> := PolynomialRing(KK, 4);

//----------------------
// Define some groups
GL4C := GL(4, KK);

// Constants
a := z^3 + z^2 + 2*z + 1;
phi := 1 + z + z^4;
// and under the action of z |-> z^2
a2 := (z^2)^3 + (z^2)^2 + 2*(z^2) + 1;
phi2 := 1 + (z^2) + (z^2)^4;

// Block diagonal matrices
g := GL4C![ phi/a,1/a,   0,0,  
            1/a,-phi/a,  0,0,
            
            0,0,         phi2/a2,1/a2,  
            0,0,         1/a2,-phi2/a2 ]; 

h := GL4C![z^2,0,  0,0,
           0,z^3,  0,0,
           
           0,0,    (z^2)^2,0,
           0,0,    0,(z^2)^3];

G := sub<GL4C | g,h>;
R := InvariantRing(G, PP);  

//----------------------
// Define some invariants
D := x0^11*x1 - 11*x0^6*x1^6 - x0*x1^11;
c4 := x0^20 + 228*x0^15*x1^5 + 494*x0^10*x1^10 - 228*x0^5*x1^15 + x1^20;
c6 := x0^30 - 522*x0^25*x1^5 - 10005*x0^20*x1^10 - 10005*x0^10*x1^20 + 522*x0^5*x1^25 + x1^30;

assert IsInvariant(D, G);
assert IsInvariant(c4, G);
assert IsInvariant(c6, G);
assert c4^3 - c6^2 eq 1728*D^5;

I44 := x0^4*y0*y1^3 - x0^3*x1*y0^4 + 3*x0^2*x1^2*y0^2*y1^2 + x0*x1^3*y1^4 - x1^4*y0^3*y1;
I17 := x0*y0^7 + 7*x0*y0^2*y1^5 - 7*x1*y0^5*y1^2 + x1*y1^7;
I71 := x0^7*y1 + 7*x0^5*x1^2*y0 + 7*x0^2*x1^5*y1 - x1^7*y0;

assert IsInvariant(I44, G);
assert IsInvariant(I71, G);
assert IsInvariant(I17, G);

w0 := I44^5;
w1 := c4*Evaluate(c4, [y0,y1,x0,x1]);
w2 := I44*(I17*I71)^2;
w3 := 12*I44^2*D*Evaluate(D, [y0,y1,x0,x1]);

//----------------------
// Define the map
A2<a,b> := AffineSpace(QQ, 2);
P3 := ProjectiveSpace(QQ, 3);

PHI := [
  b,
  -16*a^4 - 64*a^3 + 160*a^2*b - 96*a^2 + 384*a*b^2 + 80*a*b - 64*a - 256*b^2 - 224*b - 16,
  4*a^2*b - 24*a*b + 36*b,
  24*a*b - 12*b^2 + 12*b
];

PHI := map<A2 -> P3 | PHI>;
S := Image(PHI);
PHI := map<A2 -> S | DefiningPolynomials(PHI)>;

// check that the w_i satisfy the relation given by the image
assert Evaluate(DefiningPolynomial(S), [w0,w1,w2,w3]) eq 0;

//----------------------
// Define the j-maps
A<[t]> := PolynomialRing(QQ, 4);

proposed_JJ := (t[1]^2*t[2]^3)/(12*t[4]^5);

proposed_JJ_1 := -(10125000000000*t[1]^7 + 10476540000000*t[1]^6*t[2] + 
10125000000000*t[1]^6*t[3] - 5332500000000*t[1]^6*t[4] + 
1812067200*t[1]^5*t[2]^2 - 35008200000*t[1]^5*t[2]*t[3] - 
690076800000*t[1]^5*t[2]*t[4] + 19913445000000*t[1]^5*t[3]^2 - 
20295360000000*t[1]^5*t[3]*t[4] + 234832500000*t[1]^5*t[4]^2 - 
209448*t[1]^4*t[2]^3 - 42783552*t[1]^4*t[2]^2*t[3] - 61487856*t[1]^4*t[2]^2*t[4]
- 2403172800*t[1]^4*t[2]*t[3]^2 + 2229908400*t[1]^4*t[2]*t[3]*t[4] + 
17714545200*t[1]^4*t[2]*t[4]^2 - 771492600000*t[1]^4*t[3]^3 + 
149801400000*t[1]^4*t[3]^2*t[4] + 1140330600000*t[1]^4*t[3]*t[4]^2 - 
26414400000*t[1]^4*t[4]^3 + 9*t[1]^3*t[2]^4 - 1800*t[1]^3*t[2]^3*t[3] + 
3750*t[1]^3*t[2]^3*t[4] - 103680*t[1]^3*t[2]^2*t[3]^2 + 
136944*t[1]^3*t[2]^2*t[3]*t[4] + 589320*t[1]^3*t[2]^2*t[4]^2 - 
26998272*t[1]^3*t[2]*t[3]^3 + 43845840*t[1]^3*t[2]*t[3]^2*t[4] - 
897984*t[1]^3*t[2]*t[3]*t[4]^2 - 214656456*t[1]^3*t[2]*t[4]^3 + 
7116076800*t[1]^3*t[3]^4 + 123379200*t[1]^3*t[3]^3*t[4] - 
22144291200*t[1]^3*t[3]^2*t[4]^2 - 22778812800*t[1]^3*t[3]*t[4]^3 + 
1078380800*t[1]^3*t[4]^4 - 12*t[1]^2*t[2]^3*t[4]^2 + 
1296*t[1]^2*t[2]^2*t[3]^2*t[4] + 3744*t[1]^2*t[2]^2*t[3]*t[4]^2 - 
864*t[1]^2*t[2]^2*t[4]^3 - 15552*t[1]^2*t[2]*t[3]^4 - 
594000*t[1]^2*t[2]*t[3]^3*t[4] - 1139760*t[1]^2*t[2]*t[3]^2*t[4]^2 - 
597936*t[1]^2*t[2]*t[3]*t[4]^3 + 1114272*t[1]^2*t[2]*t[4]^4 + 
4940352*t[1]^2*t[3]^5 + 87891264*t[1]^2*t[3]^4*t[4] + 
174583296*t[1]^2*t[3]^3*t[4]^2 + 287122944*t[1]^2*t[3]^2*t[4]^3 + 
195658496*t[1]^2*t[3]*t[4]^4 - 17175264*t[1]^2*t[4]^5 + 
2208*t[1]*t[2]*t[3]*t[4]^4 - 1332*t[1]*t[2]*t[4]^5 - 576000*t[1]*t[3]^3*t[4]^3 -
910848*t[1]*t[3]^2*t[4]^4 - 706144*t[1]*t[3]*t[4]^5 + 109000*t[1]*t[4]^6 - 
144*t[4]^7)/(36*t[4]^5*(75000*t[1]^2 - 3*t[1]*t[2] + 600*t[1]*t[3] - 
1250*t[1]*t[4] + 4*t[4]^2));

//----------------------
// Check the moduli interpretation is what it says
JJ := (c4*Evaluate(c4, [y0,y1,x0,x1]))^3/(1728^2*(D*Evaluate(D, [y0,y1,x0,x1]))^5);
JJ_1 := (c6*Evaluate(c6, [y0,y1,x0,x1]))^2/(1728^2*(D*Evaluate(D, [y0,y1,x0,x1]))^5);

assert JJ eq Evaluate(proposed_JJ, [w0,w1,w2,w3]);
assert JJ_1 eq Evaluate(proposed_JJ_1, [w0,w1,w2,w3]);

//----------------------
// Now pull back the maps along PHI
JJ := Evaluate(proposed_JJ, DefiningPolynomials(PHI));
JJ_1 := Evaluate(proposed_JJ_1, DefiningPolynomials(PHI));
