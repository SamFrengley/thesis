/**************************************************
NUMERICAL INVARIANTS 
**************************************************/

intrinsic rho(N::RngIntElt,m::RngIntElt) -> RngIntElt
{The invariant rho(N,m) = # (x in Z/NZ : x^2 = m) }
                                                
    if GCD(N,m) ne 1 then
        return 0;
    end if;
        
    k := Valuation(N,2);
    M := N div 2^k;

    if M eq 1 then
        ret := 1;
    else
        ret := &*[1 + LegendreSymbol(m,p) : p in PrimeDivisors(M)];
    end if;
            
    if k ge 3 then
        if m mod 8 eq 1 then
            ret := 4*ret;
        else
            ret := 0;
        end if;
    elif k eq 2 then
        if m mod 4 eq 1 then
            ret := 2*ret;
        else
            ret := 0;
        end if;
    end if;

    return ret;
end intrinsic;    


/**************************************************
GENERA OF THE MODULAR CURVES X_g
**************************************************/

intrinsic XgNumbers(N::RngIntElt,r::RngIntElt) -> RngIntElt, RngIntElt, RngIntElt, RngIntElt
{When N,r are odd, computes the integers mu, e2, e3, einf}
    require (N mod 2 eq 1): "You didn't put in an odd number";
    require GCD(N,r) eq 1: "You need N,r to be coprime";
    
    if N eq 1 then
        return 1,1,1,1;
    end if;
        
    e2 := rho(N,r);
    e3 := rho(N,3*r);

    e_inf := &*[ EulerPhi(p^Valuation(N,p)) + &+[(1 + LegendreSymbol(-r, p))*EulerPhi(p^l) : l in [0..Valuation(N,p)-1] ] : p in PrimeDivisors(N) | p ne 2]; 

    if N mod 2 ne 0 then 
        mu := N^2*&*[1 + LegendreSymbol(-r,p)*1/p : p in PrimeDivisors(N) | p ne 2];
    else
        n := N div 2;
        e_inf := (EulerPhi(2^Valuation(N,2)) + &+[2*EulerPhi(2^l) : l in [0..Valuation(N,2)-1]])*e_inf;
        mu := 6*n^2*&*[1 + LegendreSymbol(-r,p)*1/p : p in PrimeDivisors(n) | p ne 2];
    end if;

    return mu, e2, e3, e_inf;
end intrinsic;


intrinsic XgPlusNumbers(N::RngIntElt,r::RngIntElt) -> FldRatElt, FldRatElt, FldRatElt, FldRatElt
{Computes the integers mu+, e2+, e3+, einf+ (these are 1/2 when N=1) }
    require N mod 2 eq 1: "Need an odd integer";
    require GCD(N,r) eq 1: "N,r need to be coprime";

    mu, e2, e3, e_inf := XgNumbers(N,r);

    return mu/2, e2/2 + ClassNumber(-4*N^2), e3/2, e_inf/2;
end intrinsic;


intrinsic GenusXg(N::RngIntElt,r::RngIntElt) -> RngIntElt
{The genus of X_g when N is odd}
    k := Valuation(N,2);
    M := N div 2^k;
    mu,e2,e3,e_inf := XgNumbers(M,r);

    require k eq 0:
      "Couldn't be bothered to impliment";
    
    if k eq 0 then 
        return [Integers()!(1 + 1/12*(mu - 3*e2 - 4*e3 - 6*e_inf))];
    end if;
end intrinsic;


intrinsic GenusXgPlus(N::RngIntElt,r::RngIntElt) -> SeqEnum[RngIntElt]
{The genus of the modular curves X_g^+ returned as a list}
    k := Valuation(N,2);
    M := N div 2^k;
    mup,e2p,e3p,e_infp := XgPlusNumbers(M,r);
    mu,e2,e3,e_inf := XgNumbers(M,r);
    
    if k eq 0 then 
        return [Integers()!(1 + 1/12*(mup - 3*e2p - 4*e3p - 6*e_infp))];
    elif k eq 1 then
        return [GenusXgPlus(M,r)[1], Integers()!(1 + 1/12*(3*mup - 3*e2p - 12*e_infp))];
    elif k eq 2 then
        if r mod 4 eq 1 then
            gw := Integers()!(1 + 1/12*(3*2^(2*k-2)*mup - 3*2*e2 - 24*e_infp));
            return [gw];
        elif r mod 4 eq 3 then
            gns := Integers()!(1 + 1/12*(2*mup - 3*ClassNumber(-4*(2*M)^2) - 4*e3 - 6*e_infp));
            gs := Integers()!(1 + 1/12*(6*mup - 3*ClassNumber(-4*(2*M)^2) - 18*e_infp));
            gw := Integers()!(1 + 1/12*(12*mup - 3*ClassNumber(-4*N^2) - 24*e_infp));
            return [gns, gs, gw];
        end if;
    elif k ge 3 then
        if r mod 8 eq 1 then
            gw := Integers()!(1 + 1/12*(3*2^(2*k-2)*mup - 3*4*e2 - 6*((2^k)*e_infp) ));
            return [gw];
        elif r mod 8 eq 3 then
            gns := Integers()!(1 + 1/12*(2^(2*k-3)*mup - 3*ClassNumber(-4*(2^(k-1)*M)^2) - 4*e3 - 6*(2^(k-2)*e_infp) ));
            gw := Integers()!(1 + 1/12*(3*2^(2*k-2)*mup - 3*ClassNumber(-4*N^2) - 6*((2^k)*e_infp) ));
            return [gns, gns, gw];
        elif r mod 8 eq 5 then
            gw := Integers()!(1 + 1/12*(3*2^(2*k-2)*mup - 6*((2^k)*e_infp) ));
            return [gw];
        elif r mod 8 eq 7 then
            gs := Integers()!(1 + 1/12*(3*2^(2*k-3)*mup - 3*ClassNumber(-4*(2^(k-1)*M)^2) - 6*(3*2^(k-2)*e_infp) ));
            gw := Integers()!(1 + 1/12*(3*2^(2*k-2)*mup - 3*ClassNumber(-4*N^2) - 6*((2^k)*e_infp) ));
            return [gs, gs, gw];
        end if;
    end if;
end intrinsic;


/**************************************************
GEOMETRIC GENUS OF WNr
**************************************************/

intrinsic GeometricGenusWNr(N::RngIntElt,r::RngIntElt) -> RngIntElt
{Computes the geometric genus of WNr as per the Theorem}
    assert GCD(N,r) eq 1;
    pg_Z := GeometricGenusZNr(N,r);

    if pg_Z eq 0 then // rational
        ret :=  0;
    end if;
    
    if N mod 2 eq 1 then 
        mu, _, _, e_inf := XgNumbers(N,r);
        ret := 1/2*(pg_Z + 1/4*(3/2*rho(N,r) + rho(N,2*r) + 2/3*rho(N,3*r) - 1/6*mu + e_inf) - 1);
	
    elif N mod 4 eq 2 then 
        M := N div 2;
        mu, _, _, e_inf := XgNumbers(M,r);
        ret := 1/2*(pg_Z + 1/4*( 2*rho(M, r) + 2/3*rho(M, 3*r) - 2/3*mu + 2*e_inf ) - 1);

    elif N mod 8 eq 4 then
	      if r mod 4 eq 1 then	    
	          M := N div 4;
	          mu, _, _, e_inf := XgNumbers(M, r);
	          ret := 1/2*(pg_Z + 1/4* ( 3*rho(M, r) + 3*e_inf - 2*mu) - 1);

        elif r mod 4 eq 3 then
            M := N div 4;
            mu, _, _, e_inf := XgNumbers(M, r);
            ret := 1/2*(pg_Z + 1/4* ( rho(M, r) + (4/3)*rho(M, 3*r) + 5*(e_inf) - 10/3*mu ) - 1);                
	      end if;

    elif N mod 8 eq 0 then
        k := [k[2] : k in Factorisation(N) | k[1] eq 2][1];
        M := N div 2^k;
        mu, _, _, e_inf := XgNumbers(M, r);
        if r mod 8 eq 1 then
            if k eq 3 then
                assert M ne 0;
                const := 6;
            else
                const := 8;
            end if;
            ret := 1/2*(pg_Z + 1/4* (const*rho(M, r) + 3*2^(k-1)*(e_inf/2) - 2^(2*k-2)*(mu/2)) - 1);

        elif r mod 8 eq 3 then
            ret := 1/2*(pg_Z + 1/4* (8/3*rho(M, 3*r) + 2^(k+1)*(e_inf/2) - 2^(2*k)/3 *(mu/2)) - 1);

        elif r mod 8 eq 5 then
            if k eq 3 then
                const := 2;
            else
                const := 0;
            end if;
            ret := 1/2*(pg_Z + 1/4* (const*rho(M, r) + 3*2^(k-1)*(e_inf/2) - 2^(2*k-2)*(mu/2)) - 1);

        elif r mod 8 eq 7 then
            ret := 1/2*(pg_Z + 1/4* (3*2^(k)*(e_inf/2) - 2^(2*k-1)*(mu/2)) - 1);
        end if;
    end if;

    return Integers()!ret;
end intrinsic;

