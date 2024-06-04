Attach("../ZNr-equations.m");
QQ := Rationals();

//////////////////////////////////////////////////
// Z(5,1) compute delta
A3<a,b,z> := AffineSpace(QQ, 3);
jj, jj_1728 := WNrModuli(5, 1 : vars:=[a,b]);

JJ := jj/1728^2;
JJ_1 := jj_1728/1728^2;
JpJ := -(JJ_1 - JJ - 1);

f := -(64*a^3 - a^2*b^2 - 152*a^2*b -48*a^2 + 2*a*b^3 + 110*a*b^2 + 166*a*b + 12*a - b^4 - 22*b^3 - 119*b^2 + 22*b - 1);

assert IsSquare(f*(JpJ^2 - 4*JJ));
S := Surface(A3, z^2 - f);

A2<s,t> := AffineSpace(QQ,2);

phi := [
  1*(2*s^2 + 15*s*t - 9*s + 12*t^2 + 11*t - 11)*(12*s^2 + 15*s*t + 11*s + 2*t^2 
    - 9*t - 11)/((2)^(2)*(s + 2*t - 1)*(2*s + t - 1)*(2*s + 3*t + 1)*(3*s + 2*t + 1)),
  -1*(2)^(2)*(s + t + 2)^(2)*(3*s^2 + 2*s + 3*t^2 + 2*t - 2)/((s + 2*t - 
     1)^(1)*(2*s + t - 1)*(2*s + 3*t + 1)*(3*s + 2*t + 1)),
  1*(s - t)*(s + t + 2)*(4*s + 7*t - 1)*(7*s + 4*t - 1)*(9*s^2 + 20*s*t + 2*s + 9*t^2
     + 2*t - 2)/((s + 2*t - 1)^(2)*(2*s + t - 1)^(2)*(2*s + 3*t + 1)*(3*s + 2*t+ 1))
];
phi := map<A2 -> S | phi>;
assert IsInvertible(phi);

jj := Evaluate(jj, DefiningPolynomials(phi));
jj_1728 := Evaluate(jj_1728, DefiningPolynomials(phi));


E1,E2 := Explode(eval Read("../Z5-r/families/5-1.m"));

//////////////////////////////////////////////////
pi := [4*t + 7*s - 1, 4*s + 7*t - 1];
aa1 := [E1[1]*pi[1]^4, E1[2]*pi[1]^6];
aa2 := [E2[1]*pi[2]^4, E2[2]*pi[2]^6];

assert forall{a : a in aa1 cat aa2 | Denominator(a) eq 1};

E1 := EllipticCurve(aa1);
E2 := EllipticCurve(aa2);
j1 := jInvariant(E1);
j2 := jInvariant(E2);

assert j1*j2 eq jj; assert jj_1728 eq (j1 - 1728)*(j2 - 1728);

DD := Discriminant(E1)*Discriminant(E2);

// Now deal with the other factors of the product of discriminants
// First compute factors of \Delta*\Delta'
PP<s,t> := PolynomialRing(Integers(),2);
DD := PP!Numerator(DD);

fcts := [fct[1] : fct in Factorisation(DD)]; S := {i : i in [1..#fcts]};
fcts := [&*[fcts[i] : i in ss] : ss in Subsets(S) | #ss ne 0];

test := [19,1129];
E1 := EllipticCurve([Evaluate(a, test) : a in aa1]);
E2 := EllipticCurve([Evaluate(a, test) : a in aa2]);

prime_tests := [p : p in PrimesInInterval(1,1000) | not p in BadPrimes(E1) cat BadPrimes(E2)];

// This loop is just a sanity check
for p in prime_tests do
    assert (TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p)) mod 5 eq 0;
end for;

// Deal with -1
assert exists{p : p in prime_tests | (TraceOfFrobenius(E1, p) - TraceOfFrobenius(QuadraticTwist(E2,-1), p)) mod 5 ne 0};

// Check that trace of Frob differ mod 5 so E and E'^{D} and E'^{-D} are not congruent
for fct in fcts do
    // Specialise to our test integers
    fct_test := Evaluate(fct, test);

    // We choose primes 1 mod 4 non-square mod 5 so that making it break at D will also work for -D. This is not required
    // it just simplifies the code.

    assert exists(q){p : p in prime_tests | (p mod 4) eq 1 and
                                            LegendreSymbol(fct_test, p) eq -1 and
                                            (TraceOfFrobenius(E1, p) mod 5) ne 0};
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, fct_test), q)) mod 5 ne 0;
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -fct_test), q)) mod 5 ne 0;
end for;

a,b := Explode(DefiningPolynomials(phi));
c6c6 := -(64*a^3 + 4488*a^2*b^2 + 2448*a^2*b - 48*a^2 + 552*a*b^4 - 14076*a*b^3 - 
4692*a*b^2 - 2484*a*b + 12*a - b^6 - 522*b^5 + 10005*b^4 + 10005*b^2 + 522*b - 1);

assert IsSquare(c6c6/(aa1[2]*aa2[2]));

//////////////////////////////////////////////////
// Z(5,2) compute delta
A3<a,b,z> := AffineSpace(QQ, 3);
jj, jj_1728 := WNrModuli(5, 2 : vars:=[a,b]);

JJ := jj/1728^2;
JJ_1 := jj_1728/1728^2;
JpJ := -(JJ_1 - JJ - 1);

f := 16*a^4 + 64*a^3 - 360*a^2*b + 96*a^2 + 216*a*b^2 - 80*a*b + 64*a + 81*b^2 + 24*b + 16;

assert IsSquare(f*(JpJ^2 - 4*JJ));
S := Surface(A3, z^2 - f);

A2<s,t> := AffineSpace(QQ,2);

phi := [
    -1*(3*s*t + 7*s + 7*t + 27)/((2)^(5)),
    1*(9*s*t - 27*s + 21*t + 1)*(9*s*t + 21*s - 27*t + 1)/((2)^(6)*(3)^(3)*(3*s + 3*t + 10)),
    1*(3*t + 7)*(s - t)*(3*s + 7)*(3*s*t + 7*s + 7*t + 123)/((2)^(8)*(3*s + 3*t + 10))
];
phi := map<A2 -> S | phi>;
assert IsInvertible(phi);

jj := Evaluate(jj, DefiningPolynomials(phi));
jj_1728 := Evaluate(jj_1728, DefiningPolynomials(phi));

E1,E2 := Explode(eval Read("../Z5-r/families/5-2.m"));

//////////////////////////////////////////////////
pi := [(3*t + 7), (3*s + 7)*(3*s*t + 7*s + 7*t + 123)];
aa1 := [E1[1]*pi[1]^4, E1[2]*pi[1]^6];
aa2 := [E2[1]*pi[2]^4, E2[2]*pi[2]^6];

assert forall{a : a in aa1 cat aa2 | Denominator(a) eq 1};

E1 := EllipticCurve(aa1);
E2 := EllipticCurve(aa2);
j1 := jInvariant(E1);
j2 := jInvariant(E2);

assert j1*j2 eq jj; assert jj_1728 eq (j1 - 1728)*(j2 - 1728);

DD := Discriminant(E1)*Discriminant(E2);

// Now deal with the other factors of the product of discriminants
// First compute factors of \Delta*\Delta'
PP<s,t> := PolynomialRing(Integers(),2);
DD := PP!Numerator(DD);

fcts := [fct[1] : fct in Factorisation(DD)]; S := {i : i in [1..#fcts]};
fcts := [&*[fcts[i] : i in ss] : ss in Subsets(S) | #ss ne 0];

test := [13,110];
E1 := EllipticCurve([Evaluate(a, test) : a in aa1]);
E2 := EllipticCurve([Evaluate(a, test) : a in aa2]);

prime_tests := [p : p in PrimesInInterval(1,200) | not p in BadPrimes(E1) cat BadPrimes(E2)];

// This loop is just a sanity check
for p in prime_tests do
    assert (TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p)) mod 5 eq 0;
end for;

// Deal with -1
assert exists{p : p in prime_tests | (TraceOfFrobenius(E1, p) - TraceOfFrobenius(QuadraticTwist(E2,-1), p)) mod 5 ne 0};

// Check that trace of Frob differ mod 5 so E and E'^{D} and E'^{-D} are not congruent
for fct in fcts do
    // Specialise to our test integers
    fct_test := Evaluate(fct, test);

    // We choose primes 1 mod 4 non-square mod 5 so that making it break at D will also work for -D. This is not required
    // it just simplifies the code.

    assert exists(q){p : p in prime_tests | (p mod 4) eq 1 and
                                            LegendreSymbol(fct_test, p) eq -1 and
                                            (TraceOfFrobenius(E1, p) mod 5) ne 0};
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, fct_test), q)) mod 5 ne 0;
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -fct_test), q)) mod 5 ne 0;
end for;

a,b := Explode(DefiningPolynomials(phi));
c6c6 := (2)*(8*a^6 + 48*a^5 - 120*a^4*b + 120*a^4 + 576*a^3*b^2 - 
300*a^3*b + 160*a^3 + 348*a^2*b^2 - 504*a^2*b + 120*a^2 - 1800*a*b^3 + 
1308*a*b^2 - 588*a*b + 48*a - 216*b^4 - 1540*b^3 - 489*b^2 - 264*b + 8);

assert IsSquare(c6c6/(aa1[2]*aa2[2]));
