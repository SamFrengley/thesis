/*
    Computational verifications for Section 3.2, when (N,r) = (3,2).
*/

QQ := Rationals();
KK<z> := CyclotomicField(3); 
PP<x0,x1,y0,y1> := PolynomialRing(KK, 4);
RR<w0,w1,w2> := PolynomialRing(QQ, 3);

//----------------------
// Define some groups
GL4 := GL(4, KK);
a := Roots(Polynomial([3,0,1]), KK)[1][1];
assert a^2 eq -3; // i.e., a = sqrt(-3)

//
g := GL4![1/a,1/(3*a),  0,0,
          6/a,-1/a,     0,0,
          
          0,0,          1/a,1/(3*a),
          0,0,          6/a,-1/a];

h := GL4![1,0,  0,0,
          0,z,  0,0,
          
          0,0,  1,0,
          0,0,  0,z^2];

// The twisted diagonal subgroup \Lambda_2
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

I22 := (18*x0*y0 + x1*y1)^2;
assert IsInvariant(I22, Lambda);

II := [ (3/2)*(3*I22^2 - c4*Evaluate(c4, [y0,y1,0,0]) - 144*D*Evaluate(D, [y0,y1,0,0])), 
        c4*Evaluate(c4, [y0,y1,0,0]),
        144*D*Evaluate(D, [y0,y1,0,0])
      ];

for inv in II do
    assert IsInvariant(inv, Lambda);
end for;

//----------------------
// Check relations
JJ := (w1/w2)^3;
JJ_1 := (w0^2 - 3*w1^2 - 3*w1*w2 - 3*w2^2)^2/(4*w2^3*(2*w0 + 3*w1 + 3*w2));

assert 1728^2*Evaluate(JJ, II) eq c4^3/D^3*Evaluate(c4^3/D^3, [y0,y1,0,0]);
assert 1728^2*Evaluate(JJ_1, II) eq (c6^2/D^3)*Evaluate(c6^2/D^3, [y0,y1,0,0]);
