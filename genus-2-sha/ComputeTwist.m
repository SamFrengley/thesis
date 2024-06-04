/**************************************************
***************************************************
****    CODE TO TWIST THE KLEIN QUARTIC        ****
***************************************************
**************************************************/

declare verbose TalkToMe, 1;


/**************************************************
Torsion subgroup polynomials, if not given in the 
APQ format
**************************************************/

function OcticPoly(f : pm:=1)
    CC<i> := ComplexField(1200);
    _<xc> := PolynomialRing(CC);
    C := HyperellipticCurve(f);
    fc := Evaluate(f, xc);

    A := AnalyticJacobian(fc);
    gg := [g : g in Generators(EndomorphismRing(A))];
    gg := [gg[2], gg[1]]; //id then RM

    M := 3*gg[1] + pm*gg[2]; //mult by a
    assert Determinant(M) eq 49;
    P := BigPeriodMatrix(A); 

    alpha := Submatrix(P*Matrix(CC,M),1,1,2,2)*Submatrix(P,1,1,2,2)^-1 ;
    alpha := Matrix(CC,2,2,[Round(Real(x)) : x in Eltseq(alpha)]); //RM should be /QQ

    p := ColumnSubmatrixRange(P, 1, 1);
    p := alpha^(-1)*p; //An RM torsion point

    pp := [FromAnalyticJacobian(n*p,A) : n in [1..3]]; //Expect to be defined over octic number field
    xpx := &+[a[1][1] + a[2][1] : a in pp];
    octic := MinimalPolynomial(xpx, 8);

    return octic;
end function;


function Poly24(f : pm:=1)
    CC<i> := ComplexField(1200);
    _<xc> := PolynomialRing(CC);
    C := HyperellipticCurve(f);
    fc := Evaluate(f, xc);

    A := AnalyticJacobian(fc);
    gg := [g : g in Generators(EndomorphismRing(A))];
    gg := [gg[2], gg[1]]; //id then RM

    M := 3*gg[1] + pm*gg[2]; //mult by a
    assert Determinant(M) eq 49;
    P := BigPeriodMatrix(A); 

    alpha := Submatrix(P*Matrix(CC,M),1,1,2,2)*Submatrix(P,1,1,2,2)^-1 ;
    alpha := Matrix(CC,2,2,[Round(Real(x)) : x in Eltseq(alpha)]); //RM should be /QQ

    p := ColumnSubmatrixRange(P, 1, 1);
    p := alpha^(-1)*p; //An RM torsion point

    pp := FromAnalyticJacobian(p,A); 
    xpx := pp[1][1] + pp[2][1];
    xx := pp[1][1] * pp[2][1];
    ply := MinimalPolynomial(xpx, 24);
    ply2 := MinimalPolynomial(xx, 24);

    return ply,ply2; 
end function;


/**************************************************
When you have more than one defining polynomial
you can get optimised repn faster
**************************************************/

function VeryOptimisedRepresentation(ff)
    K := NumberField(ff[1]);
    gg := [g*LCM([Denominator(c) : c in Coefficients(g)]) : g in ff];
    gg := [g^Matrix(Integers(), 2, 2, [1,0,0,LeadingCoefficient(g)]) div LeadingCoefficient(g) : g in gg];
    aa := &cat[[r[1] : r in Roots(g, K)] : g in gg];
    O := Order(aa : Verify:=false);
    Omax := MaximalOrder(O);
    _, Omax := OptimizedRepresentation(Omax);
    return NumberField(Omax), LLL(Omax);
end function;


/**************************************************
The action of PSL_2(F_7) on X(7).
**************************************************/

function RepresentationRho(G, GL3, z)
    assert G eq SL(2, GF(7));
    KK := BaseRing(GL3);

    S := G![0,1,-1,0]; T := G![1,1,0,1];
    G := sub<G | [S,T]>; // G but with nice generators
    a := 1 + 2*(z^3 + z^5 + z^6); 

    gg :=  GL3![[e/a : e in h] : h in [
            [z - z^6, z^2 - z^5, z^4 - z^3],
            [z^2 - z^5, z^4 - z^3, z - z^6],
            [z^4 - z^3, z - z^6, z^2 - z^5]
            ]];
    hh :=  GL3![
            [z^4, 0, 0],
            [0, z^2, 0],
            [0, 0, z]
            ];

    our_map := [<S, gg>, <T, hh>];

    return hom<G -> GL3 | our_map>;
end function;


/**************************************************
Choice of prime, chooses so that the 7th cyclotomic
polynomial is irreducible modulo p, and p does not
divide the bad primes.
**************************************************/

function ChoosePrime(badprimes, D : bigger_p:=0)
  notfound := true;
  p := 11;
  count:=0;
  while notfound do
    p := NextPrime(p);
    if (not p in badprimes) and (not (Numerator(D)*Denominator(D) mod p) eq 0) then
      _<x> := PolynomialRing(GF(p));
      if #Factorisation(x^7 - 1) eq 2 then
        if (count ge bigger_p) then
          return p;
          notfound := false;
        else
          count+:=1;
        end if;
      end if;
    end if;
  end while;
end function;


/**************************************************
Given generators P,Q for Kum(J)[3 + sqrt(2)] over
F_p chooses a symplectic basis in their span. Idea
is to take preimages on the Jacobian then comptue the
Weil pairing there.
**************************************************/

function MakeSymplectic(P, Q, pts, z_p)
    FF := BaseRing(Parent(P));
    FF := GF(#FF^2);
    Kum2 := BaseChange(Parent(P), FF);
    J2 := Jacobian(Kum2);
    P2 := Kum2!Eltseq(P);
    Q2 := Kum2!Eltseq(Q);
    PJ := Points(J2, P2)[1];
    QJ := Points(J2, Q2)[1];
    assert exists(pow){i : i in [1..6] | WeilPairing(PJ,QJ,7)^i eq z_p};
    QJ := pow*QJ;
    assert WeilPairing(PJ,QJ,7) eq z_p;
    if WeilPairing(PJ, PJ + QJ, 7) eq z_p then
        ret := PJ + QJ;
    else
        ret := PJ - QJ;
    end if;

    Pret := P;
    Qret := pow*Q;
    assert exists(PpQ){PpQ : PpQ in pts | Kum2!Eltseq(PpQ) eq Kum2!ret};

    return Pret, Qret, PpQ;
end function;


/**************************************************
Given an pp-adic approximation for a Hilbert 90 
matrix h \in GL_3, this function computes the 
twist of X(7) given by h.
**************************************************/

function ComputeTheTwist(h)
    KKp := FieldOfFractions(Parent(h[1][1]));

    P2Kp<x,y,z> := PolynomialRing(KKp, [1,1,1]);
    X7 := x^3*y + y^3*z + z^3*x;

    hh := Matrix(P2Kp, 3, 3, Eltseq(h));
    xx := hh*Matrix(P2Kp, 3, 1, [x,y,z]);
    X7_M := Evaluate(X7, Eltseq(xx));

    // Now we have a pp-adic approximation
    Qp := pAdicField(Prime(KKp), Precision(Coefficients(X7_M)[1]));
    X7_M := X7_M / Coefficients(X7_M)[1];

    // Check that in fact (as we know should be the case from the
    // twisting principle) it's actually defined over Q_p
    cc := Coefficients(X7_M); ms := Monomials(X7_M);
    cc := [Eltseq(c) : c in cc];
    for c in cc do
        assert forall{i : i in [2..#c] | Valuation(c[i]) gt 20};
    end for;

    // Coerce the equation over Q_p (now that we know this is
    // reasonable. Even better, we should be over QQ, so let's
    // spot the coefficients as rational too.
    cc := [Qp!c[1] : c in cc];
    cc := [Roots(PowerRelation(c, 1), Rationals())[1][1] : c in cc];
    d := LCM([Denominator(c) : c in cc]);
    cc := [c*d : c in cc];

    // Define the model over Q, and return.
    P2<X,Y,Z> := ProjectiveSpace(Rationals(), 2); 
    X7 := X^3*Y + Y^3*Z + Z^3*X;
    R := CoordinateRing(P2);

    ms := [Evaluate(m, [X,Y,Z]) : m in ms];
    X7_M := &+[Numerator(cc[i])*ms[i] : i in [1..#cc]];
    X7_M := Curve(P2, X7_M);

    X7_M1, mat := MinRedTernaryForm(DefiningPolynomial(X7_M));
    X7_M1 := Curve(P2, X7_M1);

    return X7_M1;
end function;



function KQT(f, g, h : bigger_p:=0)
    vprint TalkToMe, 1: "Starting twist calculation";
    QQ := Rationals();
    _<x> := PolynomialRing(QQ);
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    Kum := KummerSurface(J);

    //The Galois Action on the roots
    vprint TalkToMe, 1: "Choosing prime";
    p := ChoosePrime(BadPrimes(IntegralModel(C)), Discriminant(g) : bigger_p:=0);
    precision := 100;

    vprint TalkToMe, 1 : "Computing the Galois action";
    vtime TalkToMe, 1: GG, _, gal_data := GaloisGroup(g : Prime:=p);
    HH := NormalSubgroups(GG : IndexEqual:=6)[1]`subgroup; // K^HH = QQ(zeta)

    //p-adic roots of g ordered according to the action of GG
    rts := [GaloisRoot(g, i, gal_data : Prec:=precision, Scaled:=false) : i in [1..24]]; 
    
    KKp := Parent(rts[1]); //p-adic field we'll use

    //root of unity choice
    vprint TalkToMe, 1 : "Choosing a root of unity and computing embedding and Galois action"; 
    assert exists(z){z[1] : z in Roots(x^7 - 1, KKp) | Valuation(z[1] - 1) eq 0}; 
    z := ChangePrecision(z, precision);

    //generator for Gal(QQ(z)/QQ) = <sigma>
    assert exists(sigma){s : s in Automorphisms(KKp) | Valuation(s(z) - z^3) gt precision/2};
    
    //spotting sigma as an element of GG
    sact := [Parent(rts[1]) | 1 : r in rts];
    for r in rts do
        assert exists(i){Index(rts, t) : t in rts | Valuation(t - sigma(r)) gt precision/2};
        sact[i] := r;
    end for;

    assert exists(sigma){g : g in GG | sact eq PermuteSequence(rts, g)};    

    //The Galois action on 
    GmH, proj := quo<GG | HH>;
    phi := map<GmH -> Integers() | [<proj(sigma^i), 3^i mod 7> : i in [1..6]]>;

    function GalActmu7(tau)
        return phi(proj(tau));
    end function;


    //Choose a symplectic basis for J[3 pm sqrt2] by working mod p
    vprint TalkToMe, 1 : "Working mod p to choose a symplectic basis for K[3 pm \sqrt2]";
    Fp := GF(p); _<xp> := PolynomialRing(Fp);
    ply := MinimalPolynomial(KKp.1);
    ply := Parent(xp)![&+[u : u in [0..p-1] | Valuation(c-u) ge 1] : c in Coefficients(ply)];
    Fq := GF(p^Degree(ply)); a := Roots(ply, Fq)[1][1]; //residue field of KKp
    z_p := Evaluate(Parent(xp)![&+[u : u in [0..p-1] | Valuation(c-u) ge 1] : c in Eltseq(z)], a);

    rts_p := [Evaluate(Parent(xp)![&+[u : u in [0..p-1] | Valuation(c-u) ge 1] : c in Eltseq(r)], a) : r in rts];
    rts2_p := [r[1] : r in Roots(h, Fq)];

    Kum_p := BaseChange(Kum, Fq);
    pts_set := &join{&join{{p : p in Points(Kum_p, [1,a,b]) | 7*p eq Kum_p!0} : b in rts2_p} : a in rts_p};

    assert #pts_set eq 24;
    
    pts := [];
    for r in rts_p do 
        assert exists(pt){pt : pt in pts_set | Eltseq(pt)[2] eq r};
        Append(~pts, pt);
    end for;

    P := pts[1];
    assert exists(Q){Q : Q in pts | not Q in [n*P : n in [1..6]]};

    P,Q,PpQ := MakeSymplectic(P, Q, pts, z_p);
    pts_data := [<P,[1,0]>, <Q,[0,1]>] cat [<PseudoAddMultiple(P, Q, PpQ, n), [1,n]> : n in [1..6]];
    pts_data := &cat[[<n*pt[1], [n*pt[2][1] mod 7, n*pt[2][2] mod 7]> : n in [1..3]] : pt in pts_data];

    //And now we have a galois action on a sequence of points identified up to multiplication 
    //by -1. We have chosen P, Q (and P+Q) so that e_N(P,Q) = zeta, our chosen rt of unity.

    //GG acts by permutations of the following sequence, which is just pairs [i,j] so that a point is
    //i*P + j*Q. 
    gal_seq := [pts_data[Index([q[1] : q in pts_data], p)][2] : p in pts]; 

    P := gal_seq[Index(pts, P)];
    Q := gal_seq[Index(pts, Q)];


    //Now the cocycle calculation
    G := SL(2, GF(7));

    function GalActM(tau, pt)
        seq := PermuteSequence(gal_seq, tau);
        i := Index(seq, pt);
        return gal_seq[i];
    end function;

    function InSL2(s)
        r := GalActmu7(s);

        nP,mP := Explode(GalActM(s^(-1), P));
        nP := r*nP mod 7;

        nQ,mQ := Explode(GalActM(s^(-1), Q));
        nQ := r*nQ mod 7;

        try
            return G![nP, nQ, mP, mQ];
        catch e
            return G![nP, -nQ, mP, -mQ];
        end try;
    end function;

    vprint TalkToMe, 1 : "Computing pp-adic approximation to the cocycle";
    QQz<zeta> := CyclotomicField(7);
    GL3 := GL(3, QQz); 
    rho := RepresentationRho(G, GL3, zeta);

    function EvalMat(A)
        return Matrix(KKp, 3, 3, [Evaluate(Parent(x)!ElementToSequence(a), z) : a in Eltseq(A)] );
    end function;

    //the cocycle
    c := map<GG -> Parent(Matrix(KKp, 3, 3, [1,0,0,0,1,0,0,0,1])) | [a -> EvalMat(rho(InSL2(a))) : a in GG]>;
  
    // Choosing some nice elts of the ring of ints of the Degree 8 subfield. 
    K := NumberField(g);

    vprint TalkToMe, 1 : "First, compute degree 8 polys (can be slow, but speeds stuff up in most cases)";
    vtime TalkToMe, 1 : ply := MinimalPolynomial(&+[r[1] : r in Roots(g, K)]);
    vtime TalkToMe, 1 : ply2 := MinimalPolynomial(&+[r[1] : r in Roots(h, K)]);
    assert Degree(ply) eq 8;

    vprint TalkToMe, 1 : "Some calculations with the degree 8 subfield (LLL, optimised represenation)";
    vtime TalkToMe, 1: L,O := VeryOptimisedRepresentation([ply, ply2]);
    
    hoom := hom<L -> K | Roots(MinimalPolynomial(L.1), K)[1][1]>;
    basO := [hoom(L!o) : o in Basis(O)];
    basO := basO[2..#basO];
    basO := [Parent(x)!Eltseq(o) : o in basO];

    for count in [1..20] do
        try
            vprint TalkToMe, 1 : "computing hilbert 90";
            vtime TalkToMe, 1: h := MyHilbert90(GG, rts, c, [Index(gal_seq, P), Index(gal_seq, Q)], basO);

            KKp := FieldOfFractions(KKp);
            h := Matrix(KKp, 3, 3, Eltseq(h));
            vprint TalkToMe, 1 : "computing twists";
            tws := [ComputeTheTwist(h), ComputeTheTwist(Transpose(h)^(-1))];
            return tws;
        catch e
            assert 0 eq 0;
            vprint TalkToMe, 1 : "Hit the weird approximation failure error, just try again.";
        end try;
    end for;
end function;


intrinsic KleinQuarticTwists(f::RngUPolElt, pm::RngIntElt) -> SeqEnum 
{Given a genus 2 curve C : y^2 = f(x) with RM by sqrt 2 compute the twists of 
X(7) which parametrise elliptic curves whose 7-torsion is (symplectically or 
antisymplectically) isomorphic to Jac(C)[3 pm \sqrt 2].}
    QQ := Rationals();
    _<x> := PolynomialRing(QQ);
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    Kum := KummerSurface(J);

    require pm in {1,-1}: "pm should be a choice of sign";
    
    vprint TalkToMe, 1 : "Computing deg 24 polys";
    g,h := Poly24(f : pm:=pm);

    g := Evaluate(g, x); g := g/Coefficients(g)[25];
    h := Evaluate(h, x);

    tw1 := KQT(g, h); vprint TalkToMe, 1 : Sprint([DefiningPolynomial(c) : c in tw1]);
    return tw1;
end intrinsic;

intrinsic KleinQuarticTwists(f::RngUPolElt) -> SeqEnum
{Given a genus 2 curve C : y^2 = f(x) with RM by sqrt 2 compute the twists of 
X(7) which parametrise elliptic curves whose 7-torsion is (symplectically or 
antisymplectically) isomorphic to Jac(C)[3 \pm \sqrt 2].}
    X7_M, X7_M3 := Explode(KleinQuarticTwists(f, 1));
    X7_N, X7_N3 := Explode(KleinQuarticTwists(f, -1));

    return [X7_M, X7_M3, X7_N, X7_N3];
end intrinsic;

intrinsic KleinQuarticTwists(C::CrvHyp) -> SeqEnum
{Given a genus 2 curve C : y^2 = f(x) with RM by sqrt 2 compute the twists of 
X(7) which parametrise elliptic curves whose 7-torsion is (symplectically or 
antisymplectically) isomorphic to Jac(C)[3 \pm \sqrt 2].}
    _<x> := PolynomialRing(Rationals());
    f := Evaluate(-DefiningPolynomial(C), [x,0,1]);
    return KleinQuarticTwists(f);
end intrinsic;

intrinsic KleinQuarticTwists(APQ::SeqEnum : max_tries:=4) -> SeqEnum 
{Given a genus 2 curve C : y^2 = f(x) with RM by sqrt 2 compute the twists of 
X(7) which parametrise elliptic curves whose 7-torsion is (symplectically or 
antisymplectically) isomorphic to Jac(C)[3 pm \sqrt 2].}
    QQ := Rationals();
    _<x> := PolynomialRing(QQ);
    f := Genus2Curve(APQ);
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    Kum := KummerSurface(J);

    vprint TalkToMe, 1 : "Computing deg 24 polys. ";
    gs,hs := Explode(Poly24Alg(APQ));

    if #gs eq 0 then
      info_str := "Both repns are reducible. ";
      vprint TalkToMe, 1 : info_str;
    elif #gs eq 1 then
      info_str := "One repn is reducible. ";
      vprint TalkToMe, 1 : info_str;
    else
      info_str := "";
    end if;
    
    P<x0,x1,x2> := PolynomialRing(Rationals(), 3);

    ret := [];
    for i in [1..#gs] do
      done := false;
      counter := 0;
      while not done and counter lt max_tries do
        try
          tw := KQT(f, gs[i], hs[i] : bigger_p:=counter); vprint TalkToMe, 1 : Sprint([P!DefiningPolynomial(c) : c in tw]);
          done := true;
        catch e
          e; tw := [];
          vprintf TalkToMe, 1: "Errored out on attempt %o, try again with new prime\n", counter+1;
          counter+:=1;
        end try;
      end while;

      if not done then
        info_str := info_str cat "Ran into another error in one case. ";
      end if;
      
      ret := ret cat tw;
    end for;
    
    return ret, info_str;
end intrinsic;



/**************************************************
FUNCTIONS FOR THE MODULI INTERPRETATION OF X_M(7)
**************************************************/
function Hmap(F)
    P := Parent(F);
    dFda := Derivative(F, 1);
    dFdb := Derivative(F, 2);
    dFdc := Derivative(F, 3);
    m := [Derivative(dFda, 1), Derivative(dFda, 2), Derivative(dFda, 3),
        Derivative(dFdb, 1), Derivative(dFdb, 2), Derivative(dFdb, 3),
        Derivative(dFdc, 1), Derivative(dFdc, 2), Derivative(dFdc, 3)];
    return (-1/54)*Determinant(Matrix(P, 3, 3, m));
end function;

function c4map(F)
    P := Parent(F);
    dFda := Derivative(F, 1);
    dFdb := Derivative(F, 2);
    dFdc := Derivative(F, 3);

    H := Hmap(F);

    m := [Derivative(dFda, 1), Derivative(dFda, 2), Derivative(dFda, 3), Derivative(H, 1),
        Derivative(dFdb, 1), Derivative(dFdb, 2), Derivative(dFdb, 3), Derivative(H, 2),
        Derivative(dFdc, 1), Derivative(dFdc, 2), Derivative(dFdc, 3), Derivative(H, 3),
        Derivative(H, 1), Derivative(H, 2), Derivative(H, 3), 0];
    return (1/9)*Determinant(Matrix(P, 4, 4, m));
end function;

function c6map(F)
    P := Parent(F);
    dFda := Derivative(F, 1);
    dFdb := Derivative(F, 2);
    dFdc := Derivative(F, 3);

    H := Hmap(F);
    c4 := c4map(F);

    m := [dFda, dFdb, dFdc,
        Derivative(H, 1), Derivative(H, 2), Derivative(H, 3),
        Derivative(c4, 1), Derivative(c4, 2), Derivative(c4, 3)];
    return (1/14)*Determinant(Matrix(P, 3, 3, m));
end function;


intrinsic ModuliKleinQuarticTwist(X::CrvPln) -> MapSch
{Given a twist of the Klein quartic, outputs the j-map}
    F := DefiningPolynomial(X);
    c4 := c4map(F);
    c6 := c6map(F);

    A1 := AffineSpace(BaseRing(Parent(F)), 1);

    return map<X -> A1 | [1728*c4^3/(c4^3 - c6^2)]>;
end intrinsic;


intrinsic ModuliKleinQuarticTwist(F::RngMPolElt) -> RngMPolElt
{Given a twist of the Klein quartic, outputs the j-map}
    c4 := c4map(F);
    c6 := c6map(F);

    A1 := AffineSpace(BaseRing(Parent(F)), 1);

    return 1728*c4^3/(c4^3 - c6^2);
end intrinsic;


intrinsic CorrectQuadraticTwist(C::CrvHyp, E::CrvEll: Bound:=150) -> FldRatElt
{Given a genus 2 curve C and an elliptic curve E which are known
to have isomorphic submodules of 7-tors, up to quadratic twist this
finds a d \in QQ for which E^d and C have truely isomorphic tors}
    d := 2*3*7* &*BadPrimes(E)* &*BadPrimes(C);
    d := &*[a[1] : a in Factorisation(d)];
    dd := Divisors(d);
    dd := dd cat [-d : d in dd];

    won := false; i := 1;

    while (not won) and i le #dd do
        plus := true;

        for p in [p : p in PrimesInInterval(1,Bound) | d mod p ne 0] do
            ap_E := TraceOfFrobenius(QuadraticTwist(E,dd[i]), p);
            n1 := #Points(ChangeRing(C, GF(p)));
            n2 := #Points(ChangeRing(C, GF(p^2)));

            tr_ap_C := p + 1 - n1;
            nm_ap_C := ExactQuotient((n1^2 + n2), 2) - (p + 1)*n1 - p;
            f1 := (ap_E^2 - tr_ap_C*ap_E + nm_ap_C) mod 7;

            if f1 ne 0 then
                plus := false;
            end if;
        end for;
        won := plus;
        ret := dd[i];

        i +:= 1;
    end while;

    require won:
        "You gave stuff that wasn't congruent in the first place";

    return ret;
end intrinsic;


intrinsic CongruencesFromTwists(C::CrvHyp, XX::SeqEnum : Bound:=1000) -> SeqEnum
{Gives examples of congrunces in a quadruple of lists indexed by XX. The optional
parameter `Bound' gives the point searching bound.}
    exs := [];

    for X in XX do
        j := ModuliKleinQuarticTwist(X);
        pts := PointSearch(X, Bound);
        Es := [];
        for pt in pts do
            E := MinimalTwist(EllipticCurveWithjInvariant(Eltseq(j(pt))[1]));
            d := CorrectQuadraticTwist(C, E);
            E := QuadraticTwist(E, d);
            Append(~Es, MinimalModel(E));
        end for;

        Append(~exs, Es);
    end for;

    return exs;
end intrinsic;


intrinsic CongruencesFromTwists(APQ::SeqEnum, XX::SeqEnum : Bound:=1000) -> SeqEnum
{Gives examples of congrunces in a quadruple of lists indexed by XX. The optional
parameter `Bound' gives the point searching bound.}
    C := HyperellipticCurve(Genus2Curve(APQ));

    return CongruencesFromTwists(C, XX : Bound:=Bound);
end intrinsic;
