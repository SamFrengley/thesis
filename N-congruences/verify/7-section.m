/*
    Computational verifications for Tables 6, 7, and 8
*/
Attach("../ZNr-equations.m");

QQ := Rationals();
A3<u,v,z> := AffineSpace(QQ, 3);
P2 := ProjectiveSpace(QQ, 2);
A2 := AffineSpace(QQ, 2);

function LocallyInsolublePrimes(C)
// This has been adapted from Magma's conics functions
    aa := Coefficients(LegendrePolynomial(C));
    sol_at_infty := {Sign(aa[1]), Sign(aa[2]), Sign(aa[3])} eq {1,-1};

    bad_p := BadPrimes(C);
    if sol_at_infty eq true then
        finite_loc_insol := [];
    else 
        finite_loc_insol := [-1];
    end if;
    
    for p in bad_p do 
        if NormResidueSymbol(-aa[2]/aa[1], -aa[3]/aa[1], p) ne 1 then 
            Append(~finite_loc_insol, p);
        end if;
    end for;

    return finite_loc_insol;
end function;

function DoubleCoverOfConicToGenusOneModel(Cplus, N, r : Bound:=3)
    // From a soluble Cplus (and a bound to search for a point
    // until). Returns a degree 2 genus one model corresponding
    // to a curve on Z(12, r) mapping to Cplus on W(12, r)
    Cplus := ProjectiveClosure(Cplus);
    pt := PointSearch(Cplus, Bound)[1];
    phi := Parametrization(Cplus, pt);
    phi := DefiningPolynomials(phi);
    phi := [phi[1]/phi[3], phi[2]/phi[3]];
    f := -ZNrEquations(N, r : vars:=[phi[1], phi[2], 0])[1];
    f := Numerator(f)*Denominator(f);
    even := &*[fct[1]^(Floor(fct[2]/2)) : fct in Factorisation(f) | fct[2] gt 1];
    return Minimize(GenusOneModel(ExactQuotient(f, even^2)));
end function;

function Degree2ModelSolubleInfty(C)
    // Checks if a degree 2 model is soluble at infty
    _<t> := PolynomialRing(Rationals());
    f := Equation(C); _<X,Z> := Parent(f);
    if MonomialCoefficient(f, X^4) gt 0 or MonomialCoefficient(f, Z^4) ge 0 then
        return true;                                        // f(1,0) or f(0,1) has soln
    else
        for g in Factorisation(f) do
            gg := Evaluate(g[1], [t,1]);                     // take first patch
            K := NumberField(gg);
            r, _ := Signature(K);
            return r gt 0;                                  // if it has no real embeddings then 
                                                            // there are no points because f will always
                                                            // be negative
        end for;
    end if;
end function;

// -------------------------------------
// Check the tables in Section 5.7
// -------------------------------------

//--------------------
// (12, 5)
Z12_5 := Surface(A3, ZNrEquations(12, 5 : vars:=[u,v,z]));

table6 := eval Read("../Z12-r/curvesonZ12-r/othercurves/Z12-5.m");

for example in table6 do
    C := Curve(Z12_5, [example[1]]);
    C := Curve(IrreducibleComponents(C)[1]);
    g := Genus(C);
    irred := IsAbsolutelyIrreducible(C);

    assert g eq example[2];
    assert irred eq example[3];

    if g eq 0 and irred then
        C := Image(map<C -> P2 | Basis(-CanonicalDivisor(C))>);
        C := Conic(C); 
        assert LocallyInsolublePrimes(C) eq example[5];

    elif g eq 0 and not irred then
        assert example[4] eq "No";

    elif g eq 1 and irred then
        Cplus := Curve(A2, Evaluate(example[1], [A2.1, A2.2, 0]));
        
        if Genus(Cplus) eq 1 then
            assert example[1] eq u^3 + 3*u^2 - (v^2 + 10)*u - 3*v^2 + 6;
            //see Example 6.4
        else 
            C := DoubleCoverOfConicToGenusOneModel(Cplus, 12, 5);

            if example[4] eq "Yes" then
                assert #Points(HyperellipticCurve(C) : Bound:=10^4) ne 0;
                E := Jacobian(C); 
                assert CremonaReference(E) eq example[6];
                assert Rank(E) eq example[7];
            
            else 
                loc_insol := [];
                if not Degree2ModelSolubleInfty(C) then
                    Append(~loc_insol, -1);
                end if;
                for p in PrimeDivisors(Integers()!Discriminant(C)) do
                    if not IsLocallySoluble(C, p) then 
                        Append(~loc_insol, p);
                    end if;
                end for;

                assert loc_insol eq example[5];

                E := Jacobian(C);
                assert CremonaReference(E) eq example[6];
                assert Rank(E) eq example[7];
            end if;
        end if;

    elif g eq 1 and not irred then
        assert example[4] eq "No";
    end if;
    
end for;


//--------------------
// (12, 7)
Z12_7 := Surface(A3, ZNrEquations(12, 7 : vars:=[u,v,z]));

table7 := eval Read("../Z12-r/curvesonZ12-r/othercurves/Z12-7.m");

for example in table7 do
    C := Curve(Z12_7, [example[1]]);
    C := Curve(IrreducibleComponents(C)[1]);
    g := Genus(C);
    irred := IsAbsolutelyIrreducible(C);
    assert g eq example[2];
    assert irred eq example[3];

    if g eq 0 and irred then
        C := Image(map<C -> P2 | Basis(-CanonicalDivisor(C))>);
        C := Conic(C); 
        loc_insol := LocallyInsolublePrimes(C);
        assert loc_insol eq example[5];

    elif g eq 0 and not irred then
        assert example[4] eq "No";

    elif g eq 1 and irred then
        Cplus := Curve(A2, Evaluate(example[1], [A2.1, A2.2, 0]));
        assert Genus(Cplus) eq 0;

        C := DoubleCoverOfConicToGenusOneModel(Cplus, 12, 7);

        if example[4] eq "Yes" then
            assert #Points(HyperellipticCurve(C) : Bound:=10^4) ne 0;
            E := Jacobian(C); 
            assert CremonaReference(E) eq example[6];
            assert Rank(E) eq example[7];
        
        else 
            loc_insol := [];
            if not Degree2ModelSolubleInfty(C) then
                Append(~loc_insol, -1);
            end if;
            for p in PrimeDivisors(Integers()!Discriminant(C)) do
                if not IsLocallySoluble(C, p) then 
                    Append(~loc_insol, p);
                end if;
            end for;

            assert loc_insol eq example[5];

            E := Jacobian(C);
            assert CremonaReference(E) eq example[6];
            assert Rank(E) eq example[7];
        end if;

    elif g eq 1 and not irred then
        assert example[4] eq "No";
    end if;
    
end for;


//--------------------
// (12, 11)
Z12_11 := Surface(A3, ZNrEquations(12, 11 : vars:=[u,v,z]));

table8 := eval Read("../Z12-r/curvesonZ12-r/othercurves/Z12-11.m");

for example in table8 do
    C := Curve(Z12_11, [example[1]]);
    C := Curve(IrreducibleComponents(C)[1]);
    g := Genus(C);
    irred := IsAbsolutelyIrreducible(C);
    assert g eq example[2];
    assert irred eq example[3];
    assert irred;

    if g eq 0 then
        C := Image(map<C -> P2 | Basis(-CanonicalDivisor(C))>);
        C := Conic(C); 
        loc_insol := LocallyInsolublePrimes(C);
        assert loc_insol eq example[5];


    elif g eq 1 then
        Cplus := Curve(A2, Evaluate(example[1], [A2.1, A2.2, 0]));
        assert Genus(Cplus) eq 0;

        C := DoubleCoverOfConicToGenusOneModel(Cplus, 12, 11);

        assert example[4] eq "Yes";
        assert #Points(HyperellipticCurve(C) : Bound:=10^4) gt 0;
        E := Jacobian(C); 
        assert CremonaReference(E) eq example[6];
        assert Rank(E) eq example[7];
    end if;
    
end for;


//--------------------
// (14, 1)
Z14_1 := Surface(A3, ZNrEquations(14, 1 : vars:=[u,v,z]));

table8 := eval Read("../Z14-r/curvesonZ14-r/othercurves/Z14-1.m");

for example in table8 do
    C := Curve(Z14_1, [example[1]]);
    C := Curve(IrreducibleComponents(C)[1]);
    g := Genus(C);
    irred := IsAbsolutelyIrreducible(C);
    assert g eq example[2];
    assert irred eq example[3];
    assert irred;

    if g eq 0 then
        C := Image(map<C -> P2 | Basis(-CanonicalDivisor(C))>);
        C := Conic(C);
        loc_insol := LocallyInsolublePrimes(C);
        assert loc_insol eq example[5];

    elif g eq 1 and irred then
        Cplus := Curve(A2, Evaluate(example[1], [A2.1, A2.2, 0]));
        assert Genus(Cplus) eq 0;

        C := DoubleCoverOfConicToGenusOneModel(Cplus, 14, 1);


        assert #Points(HyperellipticCurve(C) : Bound:=10^4) ne 0;
        E := Jacobian(C);
        assert CremonaReference(E) eq example[6];
        assert Rank(E) eq example[7];
    end if;
    
end for;


//--------------------
// (14, 3)
Z14_3 := Surface(A3, ZNrEquations(14, 3 : vars:=[u,v,z]));

table8 := eval Read("../Z14-r/curvesonZ14-r/othercurves/Z14-3.m");

for example in table8 do
    C := Curve(Z14_3, [example[1]]);
    C := Curve(IrreducibleComponents(C)[1]);
    g := Genus(C);
    irred := IsAbsolutelyIrreducible(C);
    assert g eq example[2];
    assert irred eq example[3];

    if g eq 0 and irred then
        C := Image(map<C -> P2 | Basis(-CanonicalDivisor(C))>);
        C := Conic(C);
        loc_insol := LocallyInsolublePrimes(C);
        assert loc_insol eq example[5];

    elif g eq 1 and irred then
        Cplus := Curve(A2, Evaluate(example[1], [A2.1, A2.2, 0]));
        assert Genus(Cplus) eq 0;

        C := DoubleCoverOfConicToGenusOneModel(Cplus, 14, 3);

        assert #Points(HyperellipticCurve(C) : Bound:=10^4) ne 0;
        E := Jacobian(C);
        assert CremonaReference(E) eq example[6];
        assert Rank(E) eq example[7];
    end if;
    
end for;
