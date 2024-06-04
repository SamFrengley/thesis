/*
    Computational verifications for Theorem 3.3.
*/

load "../Z2_3_4/polys.m"; 
load "../Z2_3_4/iscongruent.m";
QQ := Rationals();
ZZ := Integers();

//-------------------------------
//Some helper functions

function Into2Variables(fn, vars)
    // Takes a rational function: fn \in (Q(t))(s)
    // returns a rational function: fn \in Q(t,s)
    num := Numerator(fn);
    den := Denominator(fn);
    cn,mn := CoefficientsAndMonomials(num);
    cd,md := CoefficientsAndMonomials(den);
    num := &+[Evaluate(cn[i], vars[1])*Evaluate(mn[i], [vars[2]]) : i in [1..#cn]];
    den := &+[Evaluate(cd[i], vars[1])*Evaluate(md[i], [vars[2]]) : i in [1..#cd]];
    return num/den;
end function;

// LCM of denominators of polynomial in Q(t1,...,tr)
ZZdens := func<ff | LCM([Denominator(c) : c in Coefficients(Numerator(ff))])>;

function Remove4thAnd6th(Eab)
    // Finding a pi such that pi^4 | a_4 and pi^6 | a^6
    // Input is [a_4, a_6]
    // Ouput is pi
    cand := [a[1] : a in Factorisation(Eab[1]) | a[2] ge 4];
    cand := [a : a in cand | IsDivisibleBy(Eab[2], a^6)];
    pwrs := [Max([i : i in [1..3] | IsDivisibleBy(Eab[1], a^(4*i)) and IsDivisibleBy(Eab[2], a^(6*i))]) : a in cand];
    ret := [<cand[i], pwrs[i]> : i in [1..#cand]];
    if ret eq [] then return 1; else return &*[r[1]^r[2] : r in ret]; end if;
end function;

function IsomInZZ(Eab)
    // Finding an isomorphic E/Q(t,j) such that coeffs are in Z(t,j) and
    // not too many 4th and 6th powers divide the discriminant.
    // Input: [a_4, a_6]
    // Output: new a4 and a6
    P := PolynomialRing(Integers(), 2);
    ret := Eab;
    den := LCM([SquarefreePart(Denominator(u)) : u in ret]);
    if den eq 1 then
        m := 0;
    else
        m := Min([m : m in [1..3] | not IsDivisibleBy(Denominator(ret[1]), den^(4*m)) and not IsDivisibleBy(Denominator(ret[2]), den^(6*m))]);
    end if;
    ret := [ret[1]*den^(4*m), ret[2]*den^(6*m)];
    den := LCM([ZZdens(u) : u in ret]);
    ret := [P!Numerator(ret[1]*den^4), P!Numerator(ret[2]*den^6)];
    num := Remove4thAnd6th(ret);
    return [ExactQuotient(ret[1], num^4), ExactQuotient(ret[2], num^6)];
end function;


// -------------------------------------
// Check that the maps S(N, r) \to \mathcal{W(N, r)} have the correct image
// -------------------------------------

//-------------------------------
// (N, r) = (2, 1)
P2<w0,w1,w2> := ProjectiveSpace(QQ, 2);
A4 := AffineSpace(QQ, 4);  // Think of the coordinates alpha, beta, JJ', (J-1)(J'-1)

Pa, Pb := CongPolys(2, 1);

// The surface S(2, 1) and \mathcal{W}(2, 1)
S21 := P2;
calW21 := Surface(A4, [ Evaluate(Pa, [A4.3, A4.4, A4.1, A4.2]), 
                        Evaluate(Pb, [A4.3, A4.4, A4.1, A4.2]) ]);

JJ := w1^3/(8*w0*w2^2);
JJ_1 := ((4*w0 - 3/2*w1 - w2)/w2)^2;
alpha := (4*w0 - 3/2*w1 - w2)/w2;
beta := w1/w2;

// Check that w_i really do satisfy the relation
phi := map<S21 -> calW21 | [alpha, beta, JJ, JJ_1] : Check:=true>;

//-------------------------------
// (N, r) = (4, r)
P3<w0,w1,w2,w3> := ProjectiveSpace(QQ, 3);
A5 := AffineSpace(QQ, 5);  // Think of the coordinates alpha, beta, gamma, JJ', (J-1)(J'-1)

Pa, Pb, Pg := CongPolys(4, 1);

// The surfaces S(4, r) and \mathcal{W}(4, r)
S41 := Surface(P3, 192*w0^2 - 96*w0*w1 + 128*w0*w3 + 48*w0*w2 - w2^2);
calW41 := Surface(A5, [ Evaluate(Pa, [A5.4, A5.5, A5.1, A5.2, A5.3]), 
                        Evaluate(Pb, [A5.4, A5.5, A5.1, A5.2, A5.3]),
                        Evaluate(Pg, [A5.4, A5.5, A5.1, A5.2, A5.3]) ]);

w0_2 := w0;
w1_2 := w1;
w2_2 := w3^2/w2;

JJ := w1_2^3/(8*w0_2*w2_2^2);
JJ_1 := ((4*w0_2 - 3/2*w1_2 - w2_2)/w2_2)^2;
alpha := (4*w0_2 - 3/2*w1_2 - w2_2)/w2_2;
beta := w1_2/w2_2;
gamma := (JJ*w2)/(beta*w3);

// Check that w_i really do satisfy the relation
phi := iso<S41 -> calW41 | [alpha, beta, gamma, JJ, JJ_1], 
                           [1/8*(A5.2^2/A5.3), A5.4/A5.3, A5.2*A5.3/A5.4, 1] 
                           : Check:=true, CheckInverse:=true>;


//-------------------------------
// (N, r) = (3, 1)
P3<w0,w1,w2,w3> := ProjectiveSpace(QQ, 3);
A4 := AffineSpace(QQ, 4);  // Think of the coordinates alpha, beta, JJ', (J-1)(J'-1)

Pa, Pb := CongPolys(3, 1);

// The surfaces S(3, 1) and \mathcal{W}(3, 1)
S31 := Surface(P3, 8*w0^2*w1 - 60*w0*w1^2 + 12*w0*w1*w2 + 
                    36*w0*w1*w3 - 9*w1^3 + 27*w1^2*w3 - 
                    27*w1*w3^2 + 9*w3^3);
calW31 := Surface(A4, [ Evaluate(Pa, [A4.3, A4.4, A4.1, A4.2]), 
                        Evaluate(Pb, [A4.3, A4.4, A4.1, A4.2]) ]);

JJ := w1*(4*w0^2 - 192*w0*w1 + 12*w0*w2 + 72*w0*w3 + 603*w1^2 - 
        144*w1*w2 - 918*w1*w3 + 9*w2^2 + 108*w2*w3 + 351*w3^2)/
        (36*w3^3);
JJ_1 := w1*w2^2/(4*w3^3);
alpha := (w0 + 3*w1 - 3*w3)*(8*w0*w1 - 3*w1^2 + 6*w1*w3 - 3*w3^2)/
        (12*w0*w3^2);
beta := (((3/2)*(3*w1 - w3)/w3)^2 - 3*(alpha + 1))/2;

// Check that w_i really do satisfy the relation
phi := map<S31 -> calW31 | [alpha, beta, JJ, JJ_1] : Check:=true>; 


//-------------------------------
// (N, r) = (3, 2)
P2<w0,w1,w2> := ProjectiveSpace(QQ, 2);
A4 := AffineSpace(QQ, 4);  // Think of the coordinates as alpha, beta, JJ', (J-1)(J'-1)

Pa, Pb := CongPolys(3, 2);

// The surfaces S(3, 2) and \mathcal{W}(3,2)
S32 := P2;
calW32 := Surface(A4, [ Evaluate(Pa, [A4.3, A4.4, A4.1, A4.2]), 
                        Evaluate(Pb, [A4.3, A4.4, A4.1, A4.2]) ]);

JJ := (w1/w2)^3;
JJ_1 := (w0^2 - 3*w1^2 - 3*w1*w2 - 3*w2^2)^2/
        (4*w2^3*(2*w0 + 3*w1 + 3*w2));
alpha := w1/w2;
beta := ((w0/w2)^2 - 3*(alpha^2 + alpha + 1))/2;

// Check that w_i really do satisfy the relation
phi := map<S32 -> calW32 | [alpha, beta, JJ, JJ_1] 
            : Check:=true>;


// -------------------------------------
// Check that the maps Z'(N,r) -> Z(1) are etale as claimed
// -------------------------------------
A<J1,J2,a,b,g, inv> := PolynomialRing(ZZ,6);

//-------------------------------
// N = 2,4
I := 6*J1*J2*(J1-1)*(J2-1)*(J1-J2)*inv - 1;

Pa, Pb, Pg := CongPolys(4, 1);
Pa := Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), a, b, g]);
Pb := Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), a, b, g]);
Pg := Evaluate(Pg, [J1*J2, (J1-1)*(J2-1), a, b, g]);

disc_a := Discriminant(Pa, a);
assert 1 in ideal<A | [I, disc_a]>;                             // is a unit

disc_b := Discriminant(Pb, b);
assert 1 in ideal<A | [I, Pa, disc_b]>;                         // is a unit

disc_g := Discriminant(Pg, g);
disc_g := ExactQuotient(disc_g, 2^12*3^3*J1^6*J2^6);            // scale by a unit
res := Resultant(disc_g, Pb, b);                                // element of I \cap R[a]/(P_a)
res := ExactQuotient(res, 2^6*J1^6*J2^6);                       // scale by a unit
assert res eq (J1*J2 - a^2 + 2*a - 1)^3;
assert 1 in ideal<A | [I, Pa, J1*J2 - a^2 + 2*a - 1]>;          // is a unit

// We check that the RHS of the condition on tau is a unit
// for the latter part of the proof (i.e., following Lemma 3.6) 
assert 1 in ideal<A | [I, Pa, a]>;                              // is a unit
delta := 3*(J1*J2 - (a + 1)^2);
assert 1 in ideal<A | [I, Pa, delta]>;                          // is a unit


//-------------------------------
// (N, r) = (3, 1)
xi := (J1*J2)^3 + 3*(J1*J2)^2*(J1 + J2) - 27*(J1*J2)^2 
        + 3*(J1*J2)*(J1 + J2)^2 + (J1 + J2)^3;
I := 6*J1*J2*(J1-1)*(J2-1)*(J1-J2)*xi*inv - 1;

Pa, Pb := CongPolys(3, 1);
Pa := Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), a, b]);
Pb := Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), a, b]);

disc_a := Discriminant(Pa, a);
assert 1 in ideal<A | [I, disc_a]>;                             // is a unit

disc_b := Discriminant(Pb, b);
disc_b := ExactQuotient(disc_b, 2^12*3^3*(J1-1)^6*(J2-1)^6);    // scale by a unit
res := Resultant(Pa, disc_b, a);
assert res eq -xi^2;


// We check that the RHS of the condition on tau is a unit
// for the latter part of the proof (i.e., following Lemma 3.6) 
delta_den := (b^3 - 3*(J1-1)*(J2-1)*b - 2*((J1-1)*(J2-1))^2);   // N.B. equal to (d/db) P_b
res := Resultant(delta_den, Pb, b);                             // element of I \cap R[a]/(P_a)
res := ExactQuotient(res, 2^4*3^3*(J1-1)^6*(J2-1)^6);           // scale by a unit
res := Resultant(res, Pa, a);                                   // element of I \cap R
assert res eq xi^2;                                             // is a unit

delta_num := 4*b^3 + (-10*a - 4)*b^2 - 20*(J1-1)*(J2-1)*b + 
                78*(J1-1)*(J2-1)*a + 36*J1*J2*(J1-1)*(J2-1) - 
                36*((J1-1)*(J2-1))^2 + 24*(J1-1)*(J2-1);
res := Resultant(delta_num, Pb, b);                             // element of I \cap R[a]/(P_a)
res := ExactQuotient(res, 2^8*3^2*(J1-1)^4*(J2-1)^4);           // scale by a unit
res := Resultant(res, Pa, a);                                   // element of I \cap R
res := ExactQuotient(res, (J1 - J2)^6);                         // scale by a unit
assert res eq xi^2;                                             // is a unit

//-------------------------------
// (N, r) = (3, 2)
I := 6*J1*J2*(J1-1)*(J2-1)*(J1-J2)*inv - 1;

Pa, Pb := CongPolys(3, 2);
Pa := Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), a, b]);
Pb := Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), a, b]);

disc_a := Discriminant(Pa, a);
assert 1 in ideal<A | [I, disc_a]>;                             // is a unit

disc_b := Discriminant(Pb, b);
disc_b := ExactQuotient(disc_b, 2^12*3^3*(J1-1)^6*(J2-1)^6);    // scale by a unit
assert 1 in ideal<A | [I, Pa, disc_b]>;                         // is a unit


// -------------------------------------
// Isomorphisms to X(1) \times P^1 and specialisation check
// -------------------------------------

// N.B., We confirm the claims for each pair (N,r) from the parametrisation until the
// end of the proof separately. The code behaves the same in each case.

//-------------------------------
// (N, r) = (2, 1)
// We do this one to simplify the parametrisation of Z(4,r)

//Define \mathcal{Z}(2, 1)
KJ<J1> := FunctionField(Rationals());
A3_K<J2,a,b> := AffineSpace(KJ, 3);
A1_K := AffineSpace(KJ, 1);

Pa, Pb := CongPolys(2, 1);
Pa := Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), a, b]);
Pb := Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), a, b]);

calZ21 := Curve(A3_K, [Pa,Pb]);

// Compute parametrisation
// The following parametrisation actually gives us Fisher's Hesse family
f1 := map<calZ21 -> A1_K | [(b + J1)/(a - (J1 - 1))]>;
assert IsInvertible(f1);


//-------------------------------
// (N, r) = (4, r)

// Define \mathcal{Z}(4, 1) (actually not that, a different model for it 
// using the above information about Z(2, 1) )
A2_K<t,g> := AffineSpace(KJ, 2);
A1_K<s> := AffineSpace(KJ, 1);
J2, a, b := Explode([Evaluate(p, [t]) : p in InverseDefiningPolynomials(f1)]);

Pa, Pb, Pg := CongPolys(4, 1);
assert Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), a, b, g]) eq 0;
assert Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), a, b, g]) eq 0;
Pg := Evaluate(Pg, [J1*J2, (J1-1)*(J2-1), a, b, g]);

calZ41 := Curve(A2_K, Numerator(Pg));

// Compute parametrisation
f1 := map<calZ41 -> A2_K | [t, 
            ((J1 - 1)*(J1*t^3 - 3*J1*t - 2*J1 - t^3)^2/((J1*t^2 + 2*J1*t + J1 - t^2 - 2*t)^2*J1))*g]>;
X := Image(f1);
f1 := iso<calZ41 -> X | DefiningPolynomials(f1), 
                        [t, g/((J1 - 1)*(J1*t^3 - 3*J1*t - 2*J1 - t^3)^2/((J1*t^2 + 2*J1*t + J1 - t^2 - 2*t)^2*J1))]>;

// The following parametrisation actually gives us Fisher's Hesse family
// This was found by running > X:= ProjectiveClosure(X); Parametrization(X, X![1,3,1]);
f2 := (16*J1^5*t^7 + 24*J1^5*t^6 - 66*J1^5*t^5 - 148*J1^5*t^4 - 12*J1^5*t^3 + 
        168*J1^5*t^2 + 142*J1^5*t + 36*J1^5 - 64*J1^4*t^7 - 72*J1^4*t^6 + 12*J1^4*t^5*g 
        + 198*J1^4*t^5 + 24*J1^4*t^4*g + 340*J1^4*t^4 - 24*J1^4*t^3*g + 24*J1^4*t^3 - 
        96*J1^4*t^2*g - 240*J1^4*t^2 - 84*J1^4*t*g - 166*J1^4*t - 24*J1^4*g - 36*J1^4 + 
        96*J1^3*t^7 + 72*J1^3*t^6 - 36*J1^3*t^5*g - 198*J1^3*t^5 - 58*J1^3*t^4*g - 
        236*J1^3*t^4 + 48*J1^3*t^3*g - 12*J1^3*t^3 + 150*J1^3*t^2*g + 72*J1^3*t^2 + 
        104*J1^3*t*g + 24*J1^3*t + 24*J1^3*g - 64*J1^2*t^7 - 24*J1^2*t^6 + 36*J1^2*t^5*g
        + 66*J1^2*t^5 + 44*J1^2*t^4*g + 44*J1^2*t^4 - 24*J1^2*t^3*g - J1^2*t^2*g^2 - 
        54*J1^2*t^2*g - 2*J1^2*t*g^2 - 20*J1^2*t*g - J1^2*g^2 + 16*J1*t^7 - 12*J1*t^5*g 
        - 10*J1*t^4*g + J1*t^2*g^2 - J1*t*g^3 + 2*J1*t*g^2 - J1*g^3 + 
        t*g^3)/((2)^(1)*(J1)*(J1 - 1)*(J1*t^3 - 3*J1*t - 2*J1 - t^3)*(2*J1^2*t^3 + 
        3*J1^2*t^2 - J1^2 - 4*J1*t^3 - 3*J1*t^2 + 2*t^3));
f2_inv := [(J1^2*s^4 + 6*J1^2*s^2 + 16*J1^2*s + 9*J1^2 - 2*J1*s^4 - 6*J1*s^2 - 16*J1*s + 
        s^4)/((2)^(2)*(J1 - 1)*(J1*s^3 - 3*J1*s - 2*J1 - s^3)),
        (3)*(J1)*(J1^3*s^6 - 15*J1^3*s^4 - 40*J1^3*s^3 - 45*J1^3*s^2 - 24*J1^3*s - 
        5*J1^3 - 3*J1^2*s^6 + 30*J1^2*s^4 + 80*J1^2*s^3 + 45*J1^2*s^2 + 24*J1^2*s + 
        32*J1^2 + 3*J1*s^6 - 15*J1*s^4 - 40*J1*s^3 - s^6)/((2)^(2)*(J1 - 1)*(J1*s^3 - 
        3*J1*s - 2*J1 - s^3)^(2))];
f2 := iso<X -> A1_K | [f2], f2_inv>; 

// We now check that twisting by any (\neq 1,delta) factor of \Delta_1 * \Delta_2 
// gives us non-congruent curves

// First define the curves \mathcal{E} and \mathcal{E}'
KJt<J1,t> := FunctionField(Rationals(), 2);
t_g := DefiningPolynomials(Inverse(f2)*Inverse(f1));            // equations for the parametrisation
J2 := Evaluate(J2, t_g); J2 := Into2Variables(J2, [J1, t]);     // coerce into Q(J1,t)
a := Evaluate(a, t_g); a := Into2Variables(a, [J1, t]);         // coerce into Q(J1,t)

j1 := 1728*J1; j2 := 1728*J2;
delta := 3*(J1*J2 - (a + 1)^2);
d := 3*delta*a*j1*j2/((j1 - 1728)*(j2 - 1728));
aa1 := [-27*j1/(j1 - 1728), -54*j1/(j1 - 1728)];
aa2 := [-27*d^2*j2/(j2 - 1728), -54*d^3*j2/(j2 - 1728)];
aa1 := IsomInZZ(aa1);                                           // adjust [a_4, a_6] by 4th and 6th powers resp
aa2 := IsomInZZ(aa2);                                           // adjust [a_4, a_6] by 4th and 6th powers resp

// Now with the model over Z[J1, t] factor the discriminant
PP<J1,t> := PolynomialRing(Integers(), 2);
DD := Discriminant(EllipticCurve([KJt!a : a in aa1]))*
        Discriminant(EllipticCurve([KJt!a : a in aa2]));
assert Denominator(DD) eq 1;                                    // just a check
DD := PP!Numerator(DD);                                         // coerse into Z[J1,t]

// Now we speciailise
test := [11,3];                                                 // a test pair of integers
E1 := EllipticCurve([Evaluate(a, test) : a in aa1]);            // specialisation of \mathcal{E}
E2 := EllipticCurve([Evaluate(a, test) : a in aa2]);            // specialisation of \mathcal{E}'
E1 := MinimalModel(E1);
E2 := MinimalModel(E2);
prime_tests := [p : p in PrimesInInterval(1,200) | 
                   not (p in BadPrimes(E1) cat BadPrimes(E2))]; // A set of test primes to falsify with

// This loop is just a sanity check
for p in prime_tests do
    assert (TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p)) mod 4 eq 0;
end for;

// First deal with -1 by showing trace of Frob differ mod 4 so E and E'^{-1} are not congruent
assert exists(q){p : p in prime_tests | (p mod 4) ne 1 and 
                                        TraceOfFrobenius(E1, p) mod 4 ne 0};
assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -1), q)) mod 4 ne 0;

// Now deal with the other factors of the product of discriminants
// First compute factors of \Delta*\Delta'
fcts := [fct[1] : fct in Factorisation(DD)]; S := {i : i in [1..#fcts]};
fcts := [&*[fcts[i] : i in ss] : ss in Subsets(S) | #ss ne 0]; 
// And remove the factor which is equal to delta=\Delta=\Delta' up to a square
Delta1 := Discriminant(EllipticCurve([KJt!a : a in aa1])); assert Denominator(Delta1) eq 1;
fcts := [fct : fct in fcts | not IsSquare((PP!Numerator(Delta1))*fct)];

// Check that trace of Frob differ mod 4 so E and E'^{D} and E'^{-D} are not congruent
for fct in fcts do
    // Specialise to our test integers
    fct_test := Evaluate(fct, test); 

    // We choose primes 1 mod 4 so that making it break at D will also work for -D. This is not required
    // it just simplifies the code.
    assert exists(q){p : p in prime_tests | (p mod 4) eq 1 and 
                                            LegendreSymbol(fct_test, p) eq -1 and
                                            not (TraceOfFrobenius(E1, p) mod 4) in {0,2}};
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, fct_test), q)) mod 4 ne 0;
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -fct_test), q)) mod 4 ne 0;
end for;

// Checking the claim about 44a1 and 44a2
E1 := EllipticCurve("38b2");
E2 := EllipticCurve("38b1"); 
flag, isog := IsIsogenous(E1, E2); assert flag and (Degree(isog) eq 5);
assert IsNrCongruent(4, 1, E1, E2);                             //Just checking
assert not IsNrCongruent(4, 3, E1, E2);                         //this could only fail at the tau step

// Checking the claim about 44a1 and 44a2
E1 := EllipticCurve("44a1");
E2 := EllipticCurve("44a2"); 
flag, isog := IsIsogenous(E1, E2); assert flag and (Degree(isog) eq 3);
assert IsNrCongruent(4, 3, E1, E2);                             //Just checking
assert not IsNrCongruent(4, 1, E1, E2);                         //this could only fail at the tau step


//-------------------------------
// (N, r) = (3, 1)

// Define \mathcal{Z}(3, 1)
KJ<J1> := FunctionField(Rationals());
A3_K<J2,a,b> := AffineSpace(KJ, 3);
A2_K<m,n> := AffineSpace(KJ, 2);
A1_K := AffineSpace(KJ, 1);

Pa, Pb := CongPolys(3, 1);
Pa := Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), a, b]);
Pb := Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), a, b]);

calZ31 := Curve(A3_K, [Pa,Pb]);

// Compute parametrisation
f1 := map<calZ31 -> A2_K | [-1/(J1^2 - J1*J2)*a^2 + 1/(J1 - J2)*a + 2*J2/(J1 - J2), 
                            1/2/(J1*J2 - J1 - J2 + 1)*b^2 - 3/2]>;
X := Image(f1);
f1 := map<calZ31 -> X | DefiningPolynomials(f1)>;
assert IsInvertible(f1);

// The following parametrisation actually gives us Fisher's Hesse family
// This was found by running > X:= ProjectiveClosure(X); Parametrization(X, X![1,3,1]);
f2 := (9*J1^2*m^2 - 2*J1^2*m*n + 3*J1^2*m + 2*J1^2 + 2*J1*m*n + 6*J1*m + J1*n^2 + 
        8*J1*n + 11*J1 - 4*n^2 - 8*n - 4)/((2)^(1)*(J1 - 1)*(J1*m*n - 2*J1*m - J1 + 
        2*m*n + 2*m - n^2 + 1));
f2 := map<X -> A1_K | [f2]>; 
assert IsInvertible(f2);

// We now check that twisting by any (\neq 1,delta) factor of \Delta_1 * \Delta_2 
// gives us non-congruent curves

// First define the curves \mathcal{E} and \mathcal{E}'
KJt<J1,t> := FunctionField(Rationals(), 2);
J2_a_b := DefiningPolynomials(Inverse(f2)*Inverse(f1));             // equations for the parametrisation
J2 := J2_a_b[1]; J2 := Into2Variables(J2, [J1, t]);                 // coerce into Q(J1,t)
a := J2_a_b[2]; a := Into2Variables(a, [J1, t]);                    // coerce into Q(J1,t)
b := J2_a_b[3]; b := Into2Variables(b, [J1, t]);                    // coerce into Q(J1,t)

j1 := 1728*J1; j2 := 1728*J2;
delta := -6*(2*b^3 - (5*a + 2)*b^2 - 10*(J1-1)*(J2-1)*b + 
            3*(J1-1)*(J2-1)*(13*a - 2 + 6*(J1 + J2)))/
            (b^3 - 3*(J1-1)*(J2-1)*b - 2*((J1-1)*(J2-1))^2);
d := 3*delta*j1*j2/((j1 - 1728)*(j2 - 1728));
aa1 := [-27*j1/(j1 - 1728), -54*j1/(j1 - 1728)];                    
aa2 := [-27*d^2*j2/(j2 - 1728), -54*d^3*j2/(j2 - 1728)];
aa1 := IsomInZZ(aa1);                                               // adjust [a_4, a_6] by 4th and 6th powers resp
aa2 := IsomInZZ(aa2);                                               // adjust [a_4, a_6] by 4th and 6th powers resp

// Now with the model over Z[J1, t] factor the discriminant
PP<J1,t> := PolynomialRing(Integers(), 2);
DD := Discriminant(EllipticCurve([KJt!a : a in aa1]))*
        Discriminant(EllipticCurve([KJt!a : a in aa2]));
assert Denominator(DD) eq 1;                                        // just a check
DD := PP!Numerator(DD);                                             // coerse into Z[J1,t]

// Now we speciailise
test := [11,3];                                                     // a test pair of integers
E1 := EllipticCurve([Evaluate(a, test) : a in aa1]);                // specialisation of \mathcal{E}
E2 := EllipticCurve([Evaluate(a, test) : a in aa2]);                // specialisation of \mathcal{E}'
E1 := MinimalModel(E1);
E2 := MinimalModel(E2);
prime_tests := [p : p in PrimesInInterval(1,200) | 
                   not (p in BadPrimes(E1) cat BadPrimes(E2))];     // A set of test primes to falsify with

// This loop is just a sanity check
for p in prime_tests do
    assert (TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p)) mod 3 eq 0;
end for;

// First deal with -1 by showing trace of Frob differ mod 3 so E and E'^{-1} are not congruent
assert exists(q){p : p in prime_tests | (p mod 4) ne 1 and 
                                        TraceOfFrobenius(E1, p) mod 3 ne 0};
assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -1), q)) mod 3 ne 0;

// Now deal with the other factors of the product of discriminants
// First compute factors of \Delta*\Delta'
fcts := [fct[1] : fct in Factorisation(DD)]; S := {i : i in [1..#fcts]};
fcts := [&*[fcts[i] : i in ss] : ss in Subsets(S) | #ss ne 0];

// Check that trace of Frob differ mod 3 so E and E'^{D} and E'^{-D} are not congruent
for fct in fcts do
    // Specialise to our test integers
    fct_test := Evaluate(fct, test); 

    // We choose primes 1 mod 4 and non-square mod 3 so that making it break at D 
    // will also make it break at -D. This is not required it just simplifies the code.
    assert exists(q){p : p in prime_tests | (p mod 4) eq 1 and 
                                            LegendreSymbol(fct_test, p) eq -1 and
                                            (TraceOfFrobenius(E1, p) mod 3) ne 0};
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, fct_test), q)) mod 3 ne 0;
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -fct_test), q)) mod 3 ne 0;
end for;


//-------------------------------
// (N, r) = (3, 2)

// Define \mathcal{Z}(3, 2)
KJ<J1> := FunctionField(Rationals());
A3_K<J2,a,b> := AffineSpace(KJ, 3);
A2_K<m,n> := AffineSpace(KJ, 2);
A1_K := AffineSpace(KJ, 1);

Pa, Pb := CongPolys(3, 2);
Pa := Evaluate(Pa, [J1*J2, (J1-1)*(J2-1), a, b]);
Pb := Evaluate(Pb, [J1*J2, (J1-1)*(J2-1), a, b]);

calZ32 := Curve(A3_K, [Pa,Pb]);

// Compute parametrisation
f1 := map<calZ32 -> A2_K | [a, 1/((J1-1)*(J2-1))*b^2 - 3*(a + 1)]>; 
X := Image(f1);
f1 := map<calZ32 -> X | DefiningPolynomials(f1)>;
assert IsInvertible(f1);

// The following parametrisation actually gives us Fisher's Hesse family
// This was found by running > X:= ProjectiveClosure(X); Parametrization(X, X![1,3,1]);
f2 := -(16*J1^2 + 12*J1*m^2 + 4*J1*m*n - 4*J1*m - J1*n^2 - 12*J1*n - 52*J1 - 48*m^2 
        - 16*m*n - 32*m)/((16*J1^2 + 12*J1*m^2 + 4*J1*m*n - 12*J1*m - J1*n^2 - 16*J1*n -
        68*J1 - 48*m^2 - 4*m*n - 24*m + 4*n^2 + 16*n + 16)^(1));
f2 := map<X -> A1_K | [f2]>;
assert IsInvertible(f2);

// We now check that twisting by any (\neq 1,delta) factor of \Delta_1 * \Delta_2 
// gives us non-congruent curves

// First define the curves \mathcal{E} and \mathcal{E}'
KJt<J1,t> := FunctionField(Rationals(), 2);
J2_a_b := DefiningPolynomials(Inverse(f2)*Inverse(f1));             // equations for the parametrisation
J2 := J2_a_b[1]; J2 := Into2Variables(J2, [J1, t]);                 // coerce into Q(J1,t)
b := J2_a_b[3]; b := Into2Variables(b, [J1, t]);                    // coerce into Q(J1,t)

j1 := 1728*J1; j2 := 1728*J2;
d := 3*b*j1*j2/((j1 - 1728)*(j2 - 1728));
aa1 := [-27*j1/(j1 - 1728), -54*j1/(j1 - 1728)];
aa2 := [-27*d^2*j2/(j2 - 1728), -54*d^3*j2/(j2 - 1728)];
aa1 := IsomInZZ(aa1);                                               // adjust [a_4, a_6] by 4th and 6th powers resp
aa2 := IsomInZZ(aa2);                                               // adjust [a_4, a_6] by 4th and 6th powers resp

// Now with the model over Z[J1, t] factor the discriminant
PP<J1,t> := PolynomialRing(Integers(), 2);
DD := Discriminant(EllipticCurve([KJt!a : a in aa1]))*Discriminant(EllipticCurve([KJt!a : a in aa2]));
assert Denominator(DD) eq 1;                                        // just a check
DD := PP!Numerator(DD);                                             // coerse into Z[J1,t]

// Now specialise           
test := [11,3];                                                     // a test pair of integers
E1 := EllipticCurve([Evaluate(a, test) : a in aa1]);                // specialisation of \mathcal{E}
E2 := EllipticCurve([Evaluate(a, test) : a in aa2]);                // specialisation of \mathcal{E}'
E1 := MinimalModel(E1);
E2 := MinimalModel(E2);
prime_tests := [p : p in PrimesInInterval(1,200) | 
                   not (p in BadPrimes(E1) cat BadPrimes(E2))];     // A set of test primes to falsify with

// This loop is just a sanity check
for p in prime_tests do
    assert (TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p)) mod 3 eq 0;
end for;

// First deal with -1 by showing trace of Frob differ mod 3 so E and E'^{-1} are not congruent
assert exists(q){p : p in prime_tests | (p mod 4) ne 1 and 
                                        TraceOfFrobenius(E1, p) mod 3 ne 0};
assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -1), q)) mod 3 ne 0;

// Now deal with the other factors of the product of discriminants
// First compute factors of \Delta*\Delta'
fcts := [fct[1] : fct in Factorisation(DD)]; S := {i : i in [1..#fcts]};
fcts := [&*[fcts[i] : i in ss] : ss in Subsets(S) | #ss ne 0];

// Check that trace of Frob differ mod 3 so E and E'^{D} and E'^{-D} are not congruent
for fct in fcts do
    // Specialise to our test integers
    fct_test := Evaluate(fct, test); 

    // We choose primes 1 mod 4 and non-square mod 3 so that making it break at D 
    // will also make it break at -D. This is not required it just simplifies the code.
    assert exists(q){p : p in prime_tests | (p mod 4) eq 1 and 
                                            LegendreSymbol(fct_test, p) eq -1 and
                                            (TraceOfFrobenius(E1, p) mod 3) ne 0};
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, fct_test), q)) mod 3 ne 0;
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -fct_test), q)) mod 3 ne 0;
end for;

//-------------------------------
// The final case is where (N,r)=(3,2) and beta = 0

// The parametrisation
Kt<t> := FunctionField(Rationals());
J1 := t;
J2 := 1/t;
a := 1;
b := 0;

j1 := 1728*J1; j2 := 1728*J2;
d := -2*j1*j2/((j1 - 1728)*(j2 - 1728));
aa1 := [-27*j1/(j1 - 1728), -54*j1/(j1 - 1728)];
aa2 := [-27*d^2*j2/(j2 - 1728), -54*d^3*j2/(j2 - 1728)];

// Choosing an integral model
cc1 := cInvariants(MinimalModel(EllipticCurve(aa1)))[1..2]; aa1 := [-27*cc1[1], -54*cc1[2]];
cc2 := cInvariants(MinimalModel(EllipticCurve(aa2)))[1..2]; aa2 := [-27*cc2[1], -54*cc2[2]];

// Factor the discriminant for the integral model
PP<t> := PolynomialRing(Integers());
DD := Discriminant(EllipticCurve([Kt!a : a in aa1]))*Discriminant(EllipticCurve([Kt!a : a in aa2]));
assert Denominator(DD) eq 1;                                    // Sanity check
DD := PP!Numerator(DD);                                         // Coerce

test := 11;                                                     // some test specialisation t=11
E1 := EllipticCurve([Evaluate(a, test) : a in aa1]);            // specialisation of \mathcal{E}
E2 := EllipticCurve([Evaluate(a, test) : a in aa2]);            // specialisation of \mathcal{E}
E1 := MinimalModel(E1);
E2 := MinimalModel(E2);
prime_tests := [p : p in PrimesInInterval(1,200) | not (p in BadPrimes(E1) cat BadPrimes(E2))];

// This loop is just a sanity check
for p in prime_tests do
    assert (TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p)) mod 3 eq 0;
end for;

// First deal with -1 by showing trace of Frob differ mod 3 so E and E'^{-1} are not congruent
assert exists(q){p : p in prime_tests | (p mod 4) ne 1 and 
                                        TraceOfFrobenius(E1, p) mod 3 ne 0};
assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -1), q)) mod 3 ne 0;

// Now deal with the other factors of the product of discriminants
// First compute factors of \Delta*\Delta'
fcts := [fct[1] : fct in Factorisation(DD)]; S := {i : i in [1..#fcts]};
fcts := [&*[fcts[i] : i in ss] : ss in Subsets(S) | #ss ne 0];

for fct in fcts do
    // Specialise to our test integers
    fct_test := Evaluate(fct, test); 

    // We choose primes 1 mod 4 and non-square mod 3 so that making it break at D 
    // will also make it break at -D. This is not required it just simplifies the code.
    assert exists(q){p : p in prime_tests | (p mod 4) eq 1 and 
                                            LegendreSymbol(fct_test, p) eq -1 and
                                            (TraceOfFrobenius(E1, p) mod 3) ne 0};
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, fct_test), q)) mod 3 ne 0;
    assert (TraceOfFrobenius(E1, q) - TraceOfFrobenius(QuadraticTwist(E2, -fct_test), q)) mod 3 ne 0;
end for;