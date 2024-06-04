/*
    Functions for doing elliptic Chabauty. N.B. we do not actually use these for our 
    proof of Prop 8.1 in "../verify/8-1-prop.m".

    The main function is 

    RichelotTwoCoverDescent(C : Verbose:=false)
        // Performs the methods of Section 8.1 on C, a hyperelliptic curve of genus 2 
        // which should 
        // (1) admit a richelot isogeny over Q
        // (2) be given by a model y^2 = f(x) where f(x) has degree 6 and is irreducible.
        // Returns the Q-points on C which is provably so if the second returned value is
        // true.
*/


function KS2GivingLocalPoints(g, KS2, S)
    loc_pts := [];
    for d in KS2 do
        has_pt := true; i := 1;
        while has_pt and i lt (#S + 1) do
            frakp := S[i];
            if not HasPoint(d*g, 2, frakp) then
                has_pt := false;
            end if;
            i := i + 1;
        end while;

        if has_pt then Append(~loc_pts, d); end if;
    end for;

    return loc_pts;
end function;

function LowestHeightPoint(Ed)
    // Search for points on a genus 1 curve.
    // This is not garunteed to find a point and is probably really inefficient

    found_pt := false; B := 1;
    while not found_pt do
        pts := Points(Ed : Bound:=B);
        if #pts gt 0 then
            found_pt := true;
        end if;

        if B gt 100 then
            return false, 0;
        end if;

        B := B + 10;
    end while;

    pts := [p : p in pts];
    hts := [HeightOnAmbient(p) : p in pts];
    _, i := Min(hts);

    return true, pts[i];
end function;

function LowestHeightPointsOrFail(g, d_with_loc_pts)
    // Tries to find a point on the curves
    // d*y^2 = g(x) where g has degree 4. Returns two lists
    // one of points found, and another of the d which we fail
    // to find a point.
    pts := []; fails := [];
    for d in d_with_loc_pts do
        Ed := HyperellipticCurve(d*g);
        flag, O := LowestHeightPoint(Ed);
        if flag then
            Append(~pts, <d, Coordinates(O)>);
        else
            Append(~fails, d);
        end if;
    end for;

    return pts, fails;
end function;

function ReturnAllxCoords(Ed_xmap : Verbose:=false, DegExt:=3, ClNum1:=true)
    // Takes a pair of an elliptic curve Ed/K and a map Ed -> P1. Returns the 
    // Q-points, x, on P1 such that there is a K-points in the fibre of the map
    // Ed -> P1 above x. The second value returned is false if it could not
    // be proved to be ALL the points. 
    // Uses Bruin's Elliptic Chabauty.
    all_x_cos := {};

    if ClNum1 then
        Ed, minmod := MinimalModel(Ed_xmap[1]);
        minmod := Inverse(minmod);
    else
        Ed, minmod := IntegralModel(Ed_xmap[1]);
        minmod := Inverse(minmod);
    end if;

    G, phi, _, flag := MordellWeilGroup(Ed);
    phi := (phi*minmod);

    if flag then
        rk_Ed := TorsionFreeRank(G);

        if Verbose then
            printf "The rank of the MW group here is %o \n", rk_Ed;
        end if;

        if rk_Ed lt DegExt then
            SetVerbose("EllChab", 1);
            pts, _ := Chabauty(phi, Ed_xmap[2] );

            for pt in pts do
                xx := Coordinates(Ed_xmap[2](phi(pt)));
                if xx[2] ne 0 then
                    Include(~all_x_cos, xx[1]/xx[2]);
                end if;
            end for;
            
            return all_x_cos, true;
        
        else
            printf "TROUBLE: The rank of the MW group here is too big, it's %o \n", rk_Ed;
            return {}, false;
        end if;

    else 
        printf "TROUBLE: We failed to compute the MW group here \n";
        return {}, false;
    end if;
end function;

function RichelotTwoCoverDescent(C : Verbose:=false) 
    // Performs the methods of Section 8.1 on C, a hyperelliptic curve of genus 2 
    // which should 
    // (1) admit a richelot isogeny over Q
    // (2) be given by a model y^2 = f(x) where f(x) has degree 6 and is irreducible.
    // Returns the Q-points on C which is provably so if the second returned value is
    // true.
    winning := true;
    assert Genus(C) eq 2;
    assert #RichelotIsogenousSurfaces(C) gt 0;
    _<x> := PolynomialRing(Rationals());
    
    assert IsSimplifiedModel(C) and IsIntegral(C);

    f := Evaluate(-DefiningPolynomials(C)[1], [x,0,1]);
    assert IsIrreducible(f);

    K<theta> := OptimisedRepresentation(Subfields(SplittingField(f), 3)[1][1]); //alternatively could min-red 
    O := RingOfIntegers(K);
    g := [g[1] : g in Factorisation(PolynomialRing(K)!f) | Degree(g[1]) eq 4][1];

    if Verbose then
        printf "The number field is K(theta) where theta is a root of the polynomial %o\n------------------------ \n", MinimalPolynomial(theta);
        display<X> := PolynomialRing(K);
        printf "We have the following factorisation of f(X) over the number field K \n f(X) = ( %o ) ( %o ) \n------------------------ \n", display!g, ExactQuotient(display!f, g);
    end if;


    S := Set(BadPrimes(C)) join {2}; 
    D := 1; for p in S do D := D*p; end for;
    S := [p[1] : p in Factorisation(O*D)];


    KS2, OtoKS2 := pSelmerGroup(2, Set(S));
    KS2 := [Inverse(OtoKS2)(d) : d in KS2]; // the selmer group K(S, 2)

    d_with_loc_pts := KS2GivingLocalPoints(g, KS2, S);

    if Verbose then
        printf "The d in O_{K, S}*/(O_{K, S}*)^2 such that the C_d have local points are: \n %o \n------------------------ \n", [K!d : d in d_with_loc_pts];
    end if;

    pts, fails := LowestHeightPointsOrFail(g, d_with_loc_pts);


    for d in fails do
        Ed := HyperellipticCurve(d*g);
        sel := TwoCoverDescent(Ed);
        if Verbose then
            printf "At d = %o we had to work harder and compute the fake 2-selmer group as: \n%o\n------------------------\n", K!d, sel;
        end if;

        if not (#sel eq 0) then
            printf "TROUBLE: We couldn't find a global point, nor prove there is none, at the following d: \n %o \nPossible expensive fix: adjoin a squareroot and work with a degree 6 extension.\n------------------------\n", d;
            winning := false;
        end if;
    end for;

    C_K := BaseChange(C, K);

    all_pts := Points(Scheme(C_K, [C_K.3]));

    for pp in pts do
        if Verbose then
            printf "We are doing elliptic chabauty when d = %o\n\n", K!pp[1];
        end if;
        Ed := HyperellipticCurve(pp[1]*g);
        Ed, psi := EllipticCurve(Ed, Ed!pp[2]);
        psi := [p : p in DefiningPolynomials(Inverse(psi))]; 

        psi := map<Ed -> ProjectiveSpace(Rationals(), 1) | [psi[1], psi[3]]>;
        xx, flag := ReturnAllxCoords(<Ed, psi> : Verbose:=Verbose);    
        winning := winning and flag;

        for an_x in xx do
            all_pts := all_pts join Points(Scheme(C_K, C_K.1 - an_x*C_K.3));
        end for;

        if Verbose then
            printf "------------------------\n";
        end if;
    end for;

    return all_pts, winning;
end function;