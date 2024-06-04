function CongPolys(N, r)
    assert N in [2,3,4];
    assert GCD(N, r) eq 1;

    if N eq 2 then
        _<JJ,JJ_1,a,b> := PolynomialRing(Integers(), 4);

        f1 := a^2 - JJ_1;
        f2 := b^3 - 3*JJ*b - 2*JJ*(a+1);

        return f1, f2;
    
    elif (N eq 3) and ((r mod 3) eq 1) then
        _<JJ,JJ_1,a,b> := PolynomialRing(Integers(), 4);
        JpJ := JJ - JJ_1 + 1;

        f1 := a^3 - 3*JJ*a - JJ*JpJ; 
        f2 := b^4 - 6*JJ_1*b^2 - 8*(JJ_1)^2*b - 3*(4*a + 1)*(JJ_1)^2; 

        return f1, f2;

    elif (N eq 3) and ((r mod 3) eq 2) then
        _<JJ,JJ_1,a,b> := PolynomialRing(Integers(), 4);

        f1 := a^3 - JJ;
        f2 := b^4 -6*(a + 1)*(JJ_1)*b^2 -8*(JJ_1)^2 *b - 3*(a - 1)^2*(JJ_1)^2;

        return f1, f2;
    
    elif (N eq 4) then
        _<JJ,JJ_1,a,b,g> := PolynomialRing(Integers(), 5);

        f1 := a^2 - JJ_1;
        f2 := b^3 - 3*JJ*b - 2*JJ*(a+1);
        f3 := g^4 - 6*JJ*b*g^2 - 16*JJ^2*g + 3*JJ^2*(4*JJ -  b^2);

        return f1, f2, f3;
    end if;
end function;


function TauPoly(N, r)
    if N eq 2 then
        return 0;

    elif N eq 3 and (r mod 3 eq 1) then
        _<JJ, JJ_1, a, b, t, c6_1, c6_2> := PolynomialRing(Integers(), 7);
        JpJ := JJ - JJ_1 + 1;

        delta := -6*(2*b^3 - (5*a + 2)*b^2 - 10*JJ_1*b + 3*JJ_1*(13*a - 2 + 6*JpJ))*(b^3 - 3*JJ_1*b - 2*JJ_1^2);
        return t^2 - 3*c6_1*c6_2*delta;
    
    elif N eq 3 and (r mod 3 eq 2) then 
        _<JJ, JJ_1, a, b, t, c6_1, c6_2> := PolynomialRing(Integers(), 7);
        JpJ := JJ - JJ_1 + 1;

        return t^2 - 3*c6_1*c6_2*b;
    
    elif N eq 4 and (r mod 4 eq 1) then 
        _<JJ, JJ_1, a, b, g, t, c6_1, c6_2> := PolynomialRing(Integers(), 8);

        delta := 3*(JJ - (a + 1)^2);
        return t^2 - 3*c6_1*c6_2*a*delta;
    
    elif N eq 4 and (r mod 4 eq 3) then 
        _<JJ, JJ_1, a, b, g, t, c6_1, c6_2> := PolynomialRing(Integers(), 8);

        return t^2 - 3*c6_1*c6_2*a;
    end if;
end function;