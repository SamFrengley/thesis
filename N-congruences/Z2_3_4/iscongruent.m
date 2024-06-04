function IsNrCongruent(N, r, E1, E2)
    // Uses formulae from Thm 3.3 to test if a pair of elliptic curves are (N, r)-congruent
    // where N \leq 4. Returns an error if (E1, E2) fails one of the tests in Thm 3.3.
    assert N in [2,3,4];
    assert GCD(N, r) eq 1;

    c61 := cInvariants(E1)[2];
    c62 := cInvariants(E2)[2];

    D := Discriminant(E1);

    j1 := jInvariant(E1); J1 := j1/1728;
    j2 := jInvariant(E2); J2 := j2/1728;

    JJ := J1*J2; JJ_1 := (J1 - 1)*(J2 - 1); JpJ := J1 + J2;

    assert JJ ne 0;
    assert JJ_1 ne 0;

    P<x> := PolynomialRing(Parent(j1));

    if N eq 2 then
        alphas := [a[1] : a in Roots(x^2 - JJ_1)];

        for alpha in alphas do
            betas := [a[1] : a in Roots(x^3 - 3*JJ*x - 2*JJ*(alpha+1))];
            if #betas gt 0 then
                return true;
            end if;
        end for;
    
    elif (N eq 3) and ((r mod 3) eq 1) then
        assert J1 ne J2;
        assert JJ^3 + 3*JJ^2*JpJ - 27*JJ^2 + 3*JJ*JpJ^2 + JpJ^3 ne 0;

        alphas := [a[1] : a in Roots(x^3 - 3*JJ*x - JJ*JpJ)]; 
        for alpha in alphas do
            betas := [a[1] : a in Roots( x^4 - 6*JJ_1*x^2 - 8*(JJ_1)^2*x - 3*(4*alpha + 1)*(JJ_1)^2 )]; 
            for beta in betas do
                delta := -6*(2*beta^3 - (5*alpha + 2)*beta^2 - 10*JJ_1*beta + 3*JJ_1*(13*alpha - 2 + 6*JpJ))/(beta^3 - 3*JJ_1*beta - 2*JJ_1^2);
                if #Roots(x^2 - 3*c61*c62*delta) gt 0 then
                    return true;
                end if;
            end for;
        end for;

    elif (N eq 3) and ((r mod 3) eq 2) then
        assert J1 ne J2;

        alphas := [a[1] : a in Roots(x^3 - JJ)];
        for alpha in alphas do
            betas := [a[1] : a in Roots(x^4 -6*(alpha + 1)*(JJ_1)*x^2 -8*(JJ_1)^2 *x - 3*(alpha - 1)^2*(JJ_1)^2 )];
            for beta in betas do
                if beta ne 0 and #Roots(x^2 - 3*c61*c62*beta) gt 0 then
                    return true;
                elif beta eq 0 and #Roots(x^2 + 2*c61*c62) gt 0 then
                    return true;
                end if;
            end for;
        end for;
    
    elif (N eq 4) and ((r mod 4) eq 1) then
        assert J1 ne J2;

        alphas := [a[1] : a in Roots(x^2 - JJ_1)];
        for alpha in alphas do
            betas := [a[1] : a in Roots(x^3 - 3*JJ*x - 2*JJ*(alpha+1))];
            for beta in betas do
                gammas := [a[1] : a in Roots(x^4 - 6 *JJ *beta *x^2 - 16 *JJ^2*x + 3 *JJ^2 *(4 *JJ -  beta^2) )];
                for gamma in gammas do
                    if #Roots(x^2 - 3*c61*c62*D*alpha) gt 0 then
                        return true;
                    end if;
                end for;
            end for;
        end for;

    elif (N eq 4) and ((r mod 4) eq 3) then
        assert J1 ne J2;

        alphas := [a[1] : a in Roots(x^2 - JJ_1)];
        for alpha in alphas do
            betas := [a[1] : a in Roots(x^3 - 3*JJ*x - 2*JJ*(alpha+1))];
            for beta in betas do
                gammas := [a[1] : a in Roots(x^4 - 6 *JJ *beta *x^2 - 16 *JJ^2*x + 3 *JJ^2 *(4 *JJ -  beta^2) )];
                for gamma in gammas do
                    if #Roots(x^2 - 3*c61*c62*alpha) gt 0 then
                        return true;
                    end if;
                end for;
            end for;
        end for;
    end if;

    return false;

end function;


function IsNrCongruentToTwist(N, r, E1, E2)
    // Uses formulae from Thm 3.3 to test if a pair of elliptic curves are (N, r)-congruent
    // UP TO A QUADRATIC TWIST, where N \leq 4. Returns an error if (E1, E2) fails one of 
    // the tests in Thm 3.3.
    assert N in [2,3,4];
    assert GCD(N, r) eq 1;

    c61 := cInvariants(E1)[2];
    c62 := cInvariants(E2)[2];

    D := Discriminant(E1);

    j1 := jInvariant(E1); J1 := j1/1728;
    j2 := jInvariant(E2); J2 := j2/1728;

    JJ := J1*J2; JJ_1 := (J1 - 1)*(J2 - 1); JpJ := J1 + J2;

    assert JJ ne 0;
    assert JJ_1 ne 0;

    P<x> := PolynomialRing(Parent(j1));

    twists := [];

    if N eq 2 then
        alphas := [a[1] : a in Roots(x^2 - JJ_1)];

        for alpha in alphas do
            betas := [a[1] : a in Roots(x^3 - 3*JJ*x - 2*JJ*(alpha+1))];
            if #betas gt 0 then
                return true, twists;
            end if;
        end for;
    
    elif (N eq 3) and ((r mod 3) eq 1) then
        assert J1 ne J2;
        assert JJ^3 + 3*JJ^2*JpJ - 27*JJ^2 + 3*JJ*JpJ^2 + JpJ^3 ne 0;

        alphas := [a[1] : a in Roots(x^3 - 3*JJ*x - JJ*JpJ)]; 
        for alpha in alphas do
            betas := [a[1] : a in Roots( x^4 - 6*JJ_1*x^2 - 8*(JJ_1)^2*x - 3*(4*alpha + 1)*(JJ_1)^2 )]; 
            for beta in betas do
                delta := -6*(2*beta^3 - (5*alpha + 2)*beta^2 - 10*JJ_1*beta + 3*JJ_1*(13*alpha - 2 + 6*JpJ))/(beta^3 - 3*JJ_1*beta - 2*JJ_1^2);
                t := 3*c61*c62*delta;
                if Parent(t) cmpeq Rationals() then 
                    t := Numerator(t)*Denominator(t);
                    t_nosq := Sign(t)*&*([1] cat [p[1] : p in Factorisation(t) | (p[2] mod 2) eq 1]);
                else    
                    t_nosq := t;
                end if;
                Append(~twists, t_nosq);
            end for;
        end for;

    elif (N eq 3) and ((r mod 3) eq 2) then
        assert J1 ne J2;

        alphas := [a[1] : a in Roots(x^3 - JJ)];
        for alpha in alphas do
            betas := [a[1] : a in Roots(x^4 -6*(alpha + 1)*(JJ_1)*x^2 -8*(JJ_1)^2 *x - 3*(alpha - 1)^2*(JJ_1)^2 )];
            for beta in betas do
                if beta ne 0 then
                    t := 3*c61*c62*beta;
                    if Parent(t) cmpeq Rationals() then 
                        t := Numerator(t)*Denominator(t);
                        t_nosq := Sign(t)*&*([1] cat [p[1] : p in Factorisation(t) | (p[2] mod 2) eq 1]);
                    else    
                        t_nosq := t;
                    end if;
                    Append(~twists, t_nosq);
                elif beta eq 0 then
                    t := -2*c61*c62;
                    if Parent(t) cmpeq Rationals() then 
                        t := Numerator(t)*Denominator(t);
                        t_nosq := Sign(t)*&*([1] cat [p[1] : p in Factorisation(t) | (p[2] mod 2) eq 1]);
                    else    
                        t_nosq := t;
                    end if;
                    Append(~twists, t_nosq);
                end if;
            end for;
        end for;
    
    elif (N eq 4) and ((r mod 4) eq 1) then
        assert J1 ne J2;

        alphas := [a[1] : a in Roots(x^2 - JJ_1)];
        for alpha in alphas do
            betas := [a[1] : a in Roots(x^3 - 3*JJ*x - 2*JJ*(alpha+1))];
            for beta in betas do
                gammas := [a[1] : a in Roots(x^4 - 6 *JJ *beta *x^2 - 16 *JJ^2*x + 3 *JJ^2 *(4 *JJ -  beta^2) )];
                for gamma in gammas do
                    t := 3*c61*c62*D*alpha;
                    if Parent(t) cmpeq Rationals() then 
                        t := Numerator(t)*Denominator(t);
                        t_nosq := Sign(t)*&*([1] cat [p[1] : p in Factorisation(t) | (p[2] mod 2) eq 1]);
                    else    
                        t_nosq := t;
                    end if;
                    Append(~twists, t_nosq);
                end for;
            end for;
        end for;

    elif (N eq 4) and ((r mod 4) eq 3) then
        assert J1 ne J2;

        alphas := [a[1] : a in Roots(x^2 - JJ_1)];
        for alpha in alphas do
            betas := [a[1] : a in Roots(x^3 - 3*JJ*x - 2*JJ*(alpha+1))];
            for beta in betas do
                gammas := [a[1] : a in Roots(x^4 - 6 *JJ *beta *x^2 - 16 *JJ^2*x + 3 *JJ^2 *(4 *JJ -  beta^2) )];
                for gamma in gammas do
                    t := 3*c61*c62*D*alpha;
                    if Parent(t) cmpeq Rationals() then 
                        t := Numerator(t)*Denominator(t);
                        t_nosq := Sign(t)*&*([1] cat [p[1] : p in Factorisation(t) | (p[2] mod 2) eq 1]);
                    else    
                        t_nosq := t;
                    end if;
                    Append(~twists, t_nosq);
                end for;
            end for;
        end for;
    end if;

    return #twists gt 0, twists;
end function;