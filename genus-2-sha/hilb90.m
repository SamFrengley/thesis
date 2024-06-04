intrinsic MyHilbert90(G::GrpPerm, gal::Map[GrpPerm, PowMap], c::Map[GrpPerm, GrpMat[FldFin]]) -> GrpMatElt
{Given an abstract permutation group G acting as a Galois group via the map `gal' and a 1-cochain
c: G -> GL_n(\bar K), constructs an element A \in GL_n such that c(g) = g(A) A^(-1)}

    //because I am stupid I need to reset c
    c := map<Domain(c) -> Codomain(c) | [g->c(g)^(-1) : g in Domain(c)]>;

    require G eq Domain(c) and G eq Domain(gal):
        "Need G to be the domain of the maps silly";

    barK := Domain(gal(G!1));
    p := Characteristic(barK);
    if p eq 0 then
        K := Rationals();
        gen := PrimitiveElement(barK);
        gens := [gen^i : i in [0..#G-1]];
    else
        K := GF(p);
        gens := Basis(barK, K);
    end if;

    GLn := Codomain(c);
    require BaseRing(GLn) eq barK:
        "Of course GLn must have the right base";

    d := #G;

    A := MatrixAlgebra(barK, Degree(GLn));

    hit_one := false;
    while not hit_one do
        x := [&+[Random(-80,80)*a : a in gens] : i in [1..Degree(GLn)^2]]; //some randomish matrix 
        beta := &+[ A!Eltseq(c(sigma)) * A![(gal(sigma))(a) : a in x] : sigma in G]; //form the thing in Prop X.3 Local Fields
        if Determinant(beta) ne 0 then
            hit_one := true;
            beta := GLn!Eltseq(beta);
        end if;
    end while;

    assert forall{sigma : sigma in G | GLn![(gal(sigma))(a) : a in Eltseq(beta)] eq c(sigma)^(-1)*beta};

    return beta; //sigma(beta) beta^(-1) = c_sigma
end intrinsic;


function GalAct(tau, pt, pts)
    seq := PermuteSequence(pts, tau);
    i := Index(seq, pt);
    return pts[i];
end function;

intrinsic MyHilbert90(G::GrpPerm, rts::SeqEnum[RngPadElt], c::Map[GrpPerm, AlgMat], gens::SeqEnum[RngIntElt]) -> AlgMatElt
{Given an abstract permutation group G acting as a Galois group via the map `gal' and a 1-cochain
c: G -> GL_n(\bar K), constructs an element A \in GL_n such that c(g) = g(A) A^(-1)}

    //because I am stupid I need to reset c
    c := map<Domain(c) -> Codomain(c) | [g->c(g)^(-1) : g in Domain(c)]>;

    require G eq Domain(c):
        "Need G to be the domain of the maps silly";

    barK := Parent(rts[1]);

    GLn := Codomain(c);
    d := #G;

    _<XX,YY> := PolynomialRing(barK, 2);
    mons := &cat[[XX^i * YY^j : i in [0..23]] : j in [0..20]];

    gen := rts[gens];

    hit_one := false;
    while not hit_one do
        x := [Random(0,1)*mons[Random(1,#mons)] : i in [1..Degree(GLn)^2]]; //some randomish matrix 
        beta := Matrix(barK, 3, 3, [0 : i in [1..9]]);
        for sigma in G do
            sigma_gen := [GalAct(sigma, gg, rts) : gg in gen];
            sigma_x := [Evaluate(xx, sigma_gen) : xx in x]; 

            beta +:= c(sigma) * Matrix(barK, 3, 3, sigma_x) ; //form the thing in Prop X.3 Local Fields
        end for;
        
        if Valuation(Determinant(beta)) lt Precision(rts[1])/2 then
            hit_one := true;
        end if;
    end while;

    return beta; //sigma(beta) beta^(-1) = c_sigma
end intrinsic;


intrinsic MyHilbert90(G::GrpPerm, rts::SeqEnum[RngPadElt], c::Map[GrpPerm, AlgMat], gens::SeqEnum[RngIntElt], basL::SeqEnum[RngUPolElt]) -> AlgMatElt
{Given an abstract permutation group G acting as a Galois group via the map `gal' and a 1-cochain
c: G -> GL_n(\bar K), constructs an element A \in GL_n such that c(g) = g(A) A^(-1)}

    //because I am stupid I need to reset c
    c := map<Domain(c) -> Codomain(c) | [g->c(g)^(-1) : g in Domain(c)]>;

    require G eq Domain(c):
        "Need G to be the domain of the maps silly";

    barK := Parent(rts[1]);

    GLn := Codomain(c);
    d := #G;

    PP := PolynomialRing(barK, 14);
    mons := &cat[[PP.i*PP.(j+7) : i in [1..7]] : j in [1..7]];

    gen := rts[gens];

    hit_one := false;
    while not hit_one do
        x := [Random(0,1)*mons[Random(1,#mons)] : i in [1..Degree(GLn)^2]]; //some randomish matrix 
        beta := Matrix(barK, 3, 3, [0 : i in [1..9]]);

        for sigma in G do
            sigma_gen := [GalAct(sigma, gg, rts) : gg in gen];
            basL_here := [Evaluate(f, sigma_gen[1]) : f in basL] cat [Evaluate(f, sigma_gen[2]) : f in basL];
            sigma_x := [Evaluate(xx, basL_here) : xx in x]; 

            beta +:= c(sigma) * Matrix(barK, 3, 3, sigma_x) ; //form the thing in Prop X.3 Local Fields
        end for;
        
        if Valuation(Determinant(beta)) lt Precision(rts[1])/2 then
            hit_one := true;
        end if;
    end while;

    return beta; //sigma(beta) beta^(-1) = c_sigma
end intrinsic;
