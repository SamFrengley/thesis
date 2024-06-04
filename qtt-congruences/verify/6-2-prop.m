/*
    Proving Prop 4.6.2 using 2-cover descent as described in Section 4.6. This
    code verifes the claims made.
*/

//------------
// SETUP
//------------

QQ := Rationals();
_<X> := PolynomialRing(QQ);
P1 := ProjectiveSpace(QQ, 1);

// The genus 2 curve we wish to consider
C := HyperellipticCurve(4*X^6 + 12*X^5 - 20*X^3  + 12*X + 1);
J := Jacobian(C);
assert Genus(C) eq 2;

// Our particular model for X+(ns3, ns5) and an isomorphism
// from it to C.
C2 := HyperellipticCurve(X^6 + 3*X^5 - 5*X^3 + 3*X, 1);
our_isom := map<C2 -> C | [C2.1, 2*C2.2 + (C2.3)^3, C2.3]>;
assert IsInvertible(our_isom);

// Writing out Table 4.6, the rational points on C2
table8_1 := {C2![1,1,1], C2![-2,1,1], C2![-1,0,1], C2![0,0,1], C2![1,-2,1], C2![-2,-2,1], 
             C2![1,1,0], C2![1,-1,0], C2![-1,-1,1], C2![0,-1,1], C2![3,-37,1], 
             C2![-4, -37, 1], C2![3,36,1], C2![-4,36,1]};
table8_1C := {our_isom(P) : P in table8_1};


//------------
// THE RANK
//------------
low, high := RankBounds(J);
assert low eq 2 and high eq 2;


//------------
// FACTORISATION
//------------

K<theta> := NumberField(X^3 + X^2 - 3*X + 3);
O_K := RingOfIntegers(K);
_<x> := PolynomialRing(K);

f := 4*x^6 + 12*x^5 - 20*x^3  + 12*x + 1;
f1 := 4*x^4 + 8*x^3 + (-2*theta^2 - 2*theta)*x^2 + (-2*theta^2 - 2*theta - 4)*x + theta^2 - 2*theta + 1;
f2 := x^2 + x + (theta^2 + theta - 4)/2;

assert f1*f2 eq f;

//------------
// THE SELMER GROUP
//------------

S := Set(BadPrimes(C)) join {2}; 
D := 1; for p in S do D := D*p; end for;
S := [p[1] : p in Factorisation(O_K*D)];

KS2, OtoKS2 := pSelmerGroup(2, Set(S));
KS2 := {@Inverse(OtoKS2)(d) : d in KS2@}; // Defining the selmer group as a subset of O_K
assert #KS2 eq 64;

//------------
// THE LOCAL COMPUTATIONS
//------------

// This will be the list of d \in KS2 such that dy^2 = f1(x) has a point everywhere
// locally. 
d_loc := {@ @};

for d in KS2 do
    has_loc_pt := true; i := 1;
    while has_loc_pt and i le #S do //iterate through the primes in S
        frakp := S[i];
        if not HasPoint(d*f1/4, 2, frakp) then //check if Cd' has a local point at frakp
            has_loc_pt := false;
        end if;
        i := i + 1;
    end while;

    if has_loc_pt then //if Cd' has a local point everywhere then give us back d
        Include(~d_loc, d); 
    end if;
end for;

d_from_paper := {@1, 2*theta^2 + 2*theta - 3, theta^2 + 2*theta - 1, theta^2 + theta - 3 @};
assert {d : d in d_loc} eq d_from_paper;

//------------
// DOING ELLIPTIC CHABAUTY
//------------

//------------
// When d=1
Ed := HyperellipticCurve(1*f1);
Ed, psi := EllipticCurve(Ed, Ed![0, theta-1, 1]);
psi := [p : p in DefiningPolynomials(Inverse(psi))]; 
xmap := map<Ed -> P1 | [psi[1], psi[3]]>;

Ed_min, minmod := MinimalModel(Ed);

G, phi, _, flag := MordellWeilGroup(Ed_min); 
assert flag; // This means that phi(G) is the full MW group
assert Rank(Ed_min) le 2;

pts, _ := Chabauty(phi*Inverse(minmod), xmap); // Can ignore the second entry since [E(K) : phi(G)] = 1 
xcoords1 := {xmap((phi*Inverse(minmod))(p)) : p in pts};

//------------
// When d = 2*theta^2 + 2*theta - 3
Ed := HyperellipticCurve((2*theta^2 + 2*theta - 3)*f1);
Ed, psi := EllipticCurve(Ed, Ed![-2, 3*theta + 3, 1]); 
psi := [p : p in DefiningPolynomials(Inverse(psi))]; 
xmap := map<Ed -> P1 | [psi[1], psi[3]]>;

Ed_min, minmod := IntegralModel(Ed);

flag, G, phi := PseudoMordellWeilGroup(Ed_min); 
assert flag; // This says that phi(G) has finite ODD index in G
assert Rank(Ed_min) le 2;
pts, R := Chabauty(phi*Inverse(minmod), xmap);
assert R eq 4;  // Since R is even and phi(G) has finite ODD index in G 
                // this is enough to prove the claim (see the ``Chabauty'''
                // documentation in Magma).
xcoords2 := {xmap((phi*Inverse(minmod))(p)) : p in pts};

//------------
// When d = theta^2 + 2*theta - 1
Ed := HyperellipticCurve((theta^2 + 2*theta - 1)*f1);
Ed, psi := EllipticCurve(Ed, Ed![1/2*(theta^2 + 2*theta - 1), -theta^2 + 9, 1]); 
psi := [p : p in DefiningPolynomials(Inverse(psi))]; 
xmap := map<Ed -> P1 | [psi[1], psi[3]]>;

Ed_min, minmod := MinimalModel(Ed);

G, phi, _, flag := MordellWeilGroup(Ed_min); assert flag; assert Rank(Ed_min) le 2;
pts, _ := Chabauty(phi*Inverse(minmod), xmap);
xcoords3 := {xmap((phi*Inverse(minmod))(p)) : p in pts};

//------------
// When d = theta^2 + theta - 3
Ed := HyperellipticCurve((theta^2 + theta - 3)*f1);
Ed, psi := EllipticCurve(Ed, Ed![1/2*(theta^2 - 2*theta - 1), 6*theta^2 - 9*theta + 9, 1]); 
psi := [p : p in DefiningPolynomials(Inverse(psi))]; 
xmap := map<Ed -> P1 | [psi[1], psi[3]]>;

Ed_min, minmod := MinimalModel(Ed);

G, phi, _, flag := MordellWeilGroup(Ed_min); assert flag; assert Rank(Ed_min) le 2;
pts, _ := Chabauty(phi*Inverse(minmod), xmap);
xcoords4 := {xmap((phi*Inverse(minmod))(p)) : p in pts};


//------------
// PUTTING IT ALL TOGETHER
//------------

xcoords1 := {[QQ!p[1], QQ!p[2]] : p in xcoords1};
xcoords2 := {[QQ!p[1], QQ!p[2]] : p in xcoords2};
xcoords3 := {[QQ!p[1], QQ!p[2]] : p in xcoords3};
xcoords4 := {[QQ!p[1], QQ!p[2]] : p in xcoords4};

xcoords := xcoords1 join xcoords2 join xcoords3 join xcoords4;
xcoords := [P1!xx : xx in xcoords];

all_pts := {};
xmap := map<C -> P1 | [C.1, C.3]>;

// Determine the points supported by these x-coordinates
for p in xcoords do
    all_pts := all_pts join Points(p @@ xmap);
end for;

// And finally we put it all together.
assert all_pts eq table8_1C;
