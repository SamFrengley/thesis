/**************************************************
NUMERICAL INVARIANTS
**************************************************/

function s_21(N,r)
//The number of components of E_(3,2) which become exceptional
    if N mod 4 ne 0 then
        return 1/2*ClassNumber(-4*N^2);
    else
        if r mod 4 eq 3 then
            return ClassNumber(-4*N^2);
        else
            return 0;
        end if;
    end if;
end function;  


function s_32(N,r)
//The number of components of E_(3,2) which become exceptional
    if N mod 3 ne 0 then
        return (1/2)*ClassNumber(-3*N^2);
    else
        assert r mod 3 ne 0;
        if r mod 3 eq 2 then
            return ClassNumber(-3*N^2);
        else
            return 0;
        end if;
    end if;
end function;


intrinsic X0Numbers(m::RngIntElt) -> RngIntElt, RngIntElt, RngIntElt, RngIntElt
{The elliptic points cusps etc of X_0(m)}
    mu := m*(&+[1/a : a in Divisors(m) | IsSquarefree(a)]);
    nu_inf := &+[EulerPhi(GCD(d, m div d)) : d in Divisors(m)];
    nu_3 := #{x : x in Integers(m) | x^2 + x + 1 eq 0};
    nu_2 := #{x : x in Integers(m) | x^2 + 1 eq 0};
    return Integers()!mu, nu_2, nu_3, nu_inf;
end intrinsic;


/**************************************************
SOME INTERSECTIONS ON ~ZNr
**************************************************/

intrinsic KZtil_F(N::RngIntElt,r::RngIntElt) -> SeqEnum[RngIntElt]
{Compute K.F for each F in the strict transform of Fo}
    k := Valuation(N,2);
    M := N div 2^k;

    mup,e2p,e3p,e_infp := XgPlusNumbers(M,r);
    mu,e2,e3,e_inf := XgNumbers(M,r);

    if k eq 0 then
        ret := [1/3*mup - 2*e_infp - 1/3*e3p];

    elif k eq 1 then
        ret := [1/3*mup - e_infp - 1/3*e3p, mup - 3*e_infp];

    elif k eq 2 then
        KFw := 2^(2*k-2)*mup - 3*2^(k-1)*e_infp;
        KFns := 2^(2*k-3)/3*mup - 2^(k-2)*e_infp - 2/3*e3p;
        KFs := 2^(2*k-3)*mup - 3*2^(k-2)*e_infp;
        if (r mod 4) eq 1 then
            return [KFw];
        elif (r mod 4) eq 3 then
            return [KFns, KFs, KFw];
        end if;      
        
    elif k ge 3 then
        KFw := 2^(2*k-2)*mup - 3*2^(k-1)*e_infp;
        KFns := 2^(2*k-3)/3*mup - 2^(k-2)*e_infp - 2/3*e3p;
        KFs := 2^(2*k-3)*mup - 3*2^(k-2)*e_infp;
        if (r mod 4) eq 1 then
            ret := [KFw];
        elif k ge 3 and (r mod 8) eq 3 then
            ret := [KFns, KFns, KFw];
        elif k ge 3 and (r mod 8) eq 7 then
            ret := [KFs, KFs, KFw];
        end if;
    end if;

    return [Integers()!x : x in ret];
end intrinsic;


intrinsic FgtilSquared(N::RngIntElt, r::RngIntElt) -> SeqEnum[RngIntElt]
{F.F where F is the strict transform of Fo}
    Kdot := KZtil_F(N,r);
    genera := GenusXgPlus(N,r);

    return [2*genera[i] - 2 - Kdot[i] : i in [1..#Kdot]];
end intrinsic;


/**************************************************
SOME INTERSECTIONS ON WNr
**************************************************/

intrinsic KWSquared(N::RngIntElt,r::RngIntElt) -> RngIntElt
{Self intersection of the canonical on WNr}
    k := Valuation(N,2);
    M := N div 2^k;

    if k eq 0 then
        ret := 13/2*rho(M,r) + 4*rho(M,2*r) + 2*rho(M,3*r);
    elif k eq 1 then
        ret := 9*rho(M,r) + 2*rho(M,3*r);
    elif k eq 2 and (r mod 4) eq 1 then
        ret := 14*rho(M,r);
    elif k eq 2 and (r mod 4) eq 3 then
        ret := 4*rho(M,r) + 4*rho(M,3*r);
    elif k eq 3 and (r mod 8 ) eq 1 then
        ret := 28*rho(M,r);
    elif k ge 3 and (r mod 8) eq 3 then
        ret := 8*rho(M,3*r);
    elif k eq 3 and (r mod 8) eq 5 then
        ret := 8*rho(M,r);
    elif k ge 4 and (r mod 8) eq 1 then
        ret := 36*rho(M,r);
    else
        ret := 0;
    end if;

    // at this point 2Kw^2 - Kz^2 + 2KzF - F^2 = ret
    // => Kw^2 = 1/2*Kz^2 - 3/2*KzF + sum(g - 1) + ret/2
    genera := GenusXgPlus(N,r);
    gg := &+[g - 1 : g in genera];
    
    KZ2 := KZtilSquared(N,r);
    KzF := &+KZtil_F(N,r);
    
    return Integers()!(1/2*KZ2 - 3/2*KzF + gg + ret/2);
end intrinsic;


intrinsic Cinf_KW(N::RngIntElt,r::RngIntElt) -> RngIntElt
{K.C_inf*,the intersection of C_inf* with the canonical on WNr}
    if N mod 2 eq 1 then
        ret := 1/2*rho(N,r);
    else
        ret := rho(N,r);
    end if;
    
    ret +:= Cinf_KZtil(N,r);
    ret +:= - 1/2*EulerPhi(N);
    return ret;
end intrinsic;


intrinsic KW_Fm(N::RngIntElt,r::RngIntElt,m::RngIntElt) -> FldRatElt, BoolElt
{
An upper bound on K.F_m* for each F_m on Ztil. Returns a pair x,flag where
flag is boolean which is true if F_m* is a nonsingular curve and 
K.F_m* = x
}
    mu,nu_2,nu_3,nu_inf := X0Numbers(m);

    if m mod 4 eq 3 then
        FmdotFF := ClassNumber(-4*m) + ClassNumber(-m);
    else
        FmdotFF := ClassNumber(-4*m);
    end if;
        
    ret := Integers()!(1/2*(1/3*mu - nu_inf - 1/3*nu_3 - varrho(N,m) - FmdotFF));

    if m lt N then
        flag := true;
    elif ret le -1 then
        if GeometricGenusWNr(N,r) gt 0 then
            flag := true;
        else
            ret := -1; flag := false;
        end if;
    else
        flag := false;
    end if;
    
    return ret, flag;
end intrinsic;


intrinsic KW_Fg(N::RngIntElt,r::RngIntElt) -> SeqEnum[RngIntElt]
{K.Fg* for each g}
    k := Valuation(N,2);
    M := N div 2^k;
        
    Ktil := KZtil_F(N,r);
    genera := GenusXgPlus(N,r);

    // the following are giving KZtil.Fg_til - KZo.Fg_o
    if k eq 0 then
        ret := [3/2*rho(M,r) + rho(M,2*r) + 1/2*rho(M,3*r)];
    elif k eq 1 then
        ret := [1/2*rho(M,3*r) + 3/2*rho(M,r), 1/2*rho(M,r)];
    elif k eq 2 and (r mod 4) eq 1 then
        ret := [3*rho(M,r)];
    elif k eq 2 and (r mod 4) eq 3 then
        ret := [rho(M,3*r), rho(M,r)];
    elif k eq 3 and (r mod 8 ) eq 1 then
        ret := [6*rho(M,r)];
    elif k ge 3 and (r mod 8) eq 3 then
        ret := [rho(M,3*r),rho(M,3*r),0];
    elif k eq 3 and (r mod 8) eq 5 then
        ret := [2*rho(M,r)];
    elif k ge 4 and (r mod 8) eq 1 then
        ret := [8*rho(M,r)];
    else
        ret := [0 : i in [1..#Ktil]];
    end if;

    // By projection we have KW.Fg* = 1/2*(KZo - FFo).(Fg_o) = 1/2*(KZo.Fgo - Fgo^2)
    // By adjunction this is KZo.Fgo - (g - 1)
    return [2*((Ktil[i] - ret[i]) - (genera[i] - 1)) : i in [1..#ret]];    
end intrinsic;


/**************************************************
Intersection numbers on _WNr
**************************************************/
intrinsic barCinf_barKW(N::RngIntElt,r::RngIntElt) -> RngIntElt
{K .\bar(C)_inf,the intersection of \bar(C)_inf with the canonical on _WNr}
    return Cinf_KW(N,r) -1/2*(&+[rho(d,-r)*EulerPhi(N div d) : d in Divisors(N) | d ne 1]);
end intrinsic;

/**************************************************
Intersection numbers on WNro
**************************************************/

function rho_2(N,r)
    if N mod 4 eq 2 then
        return rho(N div 2, r);
    else
        return 0;
    end if;
end function;

function rho_3(N,r)
    if (N mod 9 eq 3) or (N mod 9 eq 6) then
        return rho(N div 3, r);
    else
        return 0;
    end if;
end function;

intrinsic KWoSquared(N::RngIntElt,r::RngIntElt) -> RngIntElt
{
Self intersection of canonical on WNr^circ
}
    require GeometricGenusWNr(N,r) ne 0: "Your surface is rational.";
    
    KW2 := KWSquared(N,r);
    sum_term := &+([0] cat [Floor(d/2)*rho(d,-r)*EulerPhi(N div d)/2 : d in Divisors(N) | d ne 1]);
    // s1 := 1/4* &+([0] cat [(d - 1)*rho(d,-r)*EulerPhi(N div d) : d in Divisors(N) | d ne 1 and (d mod 2) eq 1]);
    // s2 := 1/4* &+([0] cat [d*rho(d,-r)*EulerPhi(N div d) : d in Divisors(N) | d ne 1 and (d mod 2) eq 0]);

    g0X0 := {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,29,31,32,35,36,39,41,47,49,50,59,71};
    g0X0 := {m : m in g0X0 | GCD(m,N) eq 1};
    exceptional_X0 := {m : m in g0X0 | KW_Fm(N,r,m) eq -1};

    // The number of exceptional F_m,lamdba
    frak_m := 1/2* &+([0] cat [rho(N,r*m) : m in exceptional_X0]);
    
    // Add the number of components of E_(2,1)* meeting an exceptional X0
    // This will just be the number of X0(5), but good to know anyway
    // s3 +:= &+([0] cat [1/2*rho(N,r*m)*NumberPrimitiveRepsx2y2(m) : m in exceptional_X0]);

    ret := KW2 + sum_term + s_21(N,r) + s_32(N,r) + frak_m + 1/2*rho(N,3*r) + 1/2*rho(N,5*r) + + 1/2*rho_2(N,r) + rho_2(N,2*r) + 3/2*rho_3(N,r);
    return Integers()!ret;
end intrinsic;

