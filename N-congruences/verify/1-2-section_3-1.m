/*
    Computational verifications for Section 3.2, when (N,r) = (3,1).
*/

QQ := Rationals();
KK<z> := CyclotomicField(3); 
PP<x0,x1,y0,y1> := PolynomialRing(KK, 4);
RR<w0,w1,w2,w3> := PolynomialRing(QQ, 4);

//----------------------
// Define some groups
GL4 := GL(4, KK); 
a := Roots(Polynomial([3,0,1]), KK)[1][1];
assert a^2 eq -3; // i.e., a = sqrt(-3)

// Define block diagonal matrices
g := GL4![1/a,(1/3)/a,  0,0,
          6/a,-1/a,     0,0,
          
          0,0,          1/a,(1/3)/a,
          0,0,          6/a,-1/a];

h := GL4![1,0,  0,0,
          0,z,  0,0,
          
          0,0,  1,0,
          0,0,  0,z];

// The diagonal subgroup \Lambda_1
Lambda := sub<GL4 | g,h>;
Lambda_tilde := CommutatorSubgroup(Lambda);

//----------------------
// Define some invariants

D := -(27*x0^4 + x0*x1^3);
c4 := -(216*x0^3*x1 - x1^4);
c6 := 5832*x0^6 - 540*x0^3*x1^3 - x1^6;

assert IsInvariant(D, Lambda_tilde);
assert IsInvariant(c4, Lambda_tilde);
assert IsInvariant(c6, Lambda_tilde);
assert c4^3 - c6^2 eq 1728*D^3;

I22 := 3*(54*x0^2*y0^2 + x0*x1*y1^2 + x1^2*y0*y1);
I66 := 3^6*(x0*y1 - x1*y0)^6;

assert IsInvariant(I22, Lambda);
assert IsInvariant(I66, Lambda);

II := [ 12*I66, 
        4*I22^3,
        c6*Evaluate(c6, [y0,y1,0,0]),
        144*I22*D*Evaluate(D, [y0,y1,0,0])
      ];

for inv in II do
    assert IsInvariant(inv, Lambda);
end for;

//----------------------
// Check relations
f := 8*w0^2*w1 - 60*w0*w1^2 + 12*w0*w1*w2 + 36*w0*w1*w3 - 9*w1^3 + 27*w1^2*w3 - 27*w1*w3^2 + 9*w3^3;
JJ := w1*(4*w0^2 - 192*w0*w1 + 12*w0*w2 + 72*w0*w3 + 603*w1^2 - 144*w1*w2 - 918*w1*w3 + 9*w2^2 + 108*w2*w3 + 351*w3^2)/(36*w3^3);
JJ_1 := w1*w2^2/(4*w3^3);

assert Evaluate(f, II) eq 0;
assert 1728^2*Evaluate(JJ, II) eq c4^3/D^3*Evaluate(c4^3/D^3, [y0,y1,0,0]);
assert 1728^2*Evaluate(JJ_1, II) eq (c6^2/D^3)*Evaluate(c6^2/D^3, [y0,y1,0,0]);
