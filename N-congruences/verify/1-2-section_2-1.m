/*
    Computational verifications for Section 3.2, when (N,r) = (2,1).
*/
KK<z> := CyclotomicField(4);
PP<x0,x1,y0,y1> := PolynomialRing(KK, 4);
RR<w0,w1,w2> := PolynomialRing(KK, 3);
//----------------------
// Define some groups
GL4 := GL(4, KK);

// Block diagonal matrices
g := [1/2,1/16,  0,0,
      12,-1/2,   0,0,
      
      0,0,       1/2,1/16,
      0,0,       12,-1/2];

h := [1,0,  0,0,
      0,-1, 0,0,
      
      0,0,  1,0,
      0,0,  0,-1];

g := GL4![x : x in g]; //scale by zeta_4
h := GL4![x : x in h]; //scale by zeta_4

// The diagonal subgroup \Lambda_1
Lambda := sub<GL4 | h,g>;
Lambda_tilde := CommutatorSubgroup(Lambda);

//----------------------
// Define some invariants
D := x0*(64*x0^2 - x1^2);
c4 := 192*x0^2 + x1^2;
c6 := x1*(576*x0^2 - x1^2);

assert IsInvariant(D, Lambda_tilde);
assert IsInvariant(c4, Lambda_tilde);
assert IsInvariant(c6, Lambda_tilde);
assert c4^3 - c6^2 eq 1728*D^2;

I11 := 192*x0*y0 + x1*y1;

assert IsInvariant(I11, Lambda);

II := [ I11^3, 
        2*I11*c4*Evaluate(c4, [y0,y1,0,0]),
        1728*D*Evaluate(D, [y0,y1,0,0])
      ];

for inv in II do
    assert IsInvariant(inv, Lambda);
end for;

//----------------------
// Check relations
JJ := w1^3/(8*w0*w2^2);
JJ_1 := ((4*w0 - 3/2*w1 - w2)/w2)^2;

assert 1728^2*Evaluate(JJ, II) eq c4^3/D^2*Evaluate(c4^3/D^2, [y0,y1,0,0]);
assert 1728^2*Evaluate(JJ_1, II) eq (c6^2/D^2)*Evaluate(c6^2/D^2, [y0,y1,0,0]);
