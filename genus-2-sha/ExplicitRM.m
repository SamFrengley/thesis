declare verbose TalkToMe, 1;

intrinsic Genus2Curve(APQ::SeqEnum) -> SeqEnum, RngUPolElt
{Theorem 5.1 (Bending)}
    if Type(Parent(APQ[1])) cmpeq RngInt then
        A,P,Q := Explode([Rationals()!a : a in APQ]);
    else
        A,P,Q := Explode(APQ);
    end if;

    PQ<x> := PolynomialRing(Parent(A));
    R := 4*P;
    B := (Q*(P*A-Q)+4*P^2+1)/P^2;
    C := 4*(P*A-Q)/P;
    f := x^6 + A*x^5 + (A^2*P^3 - 3*A*P^2*Q + A*P*Q + 4*P^3 
        + 4*P^2 + 2*P*Q^2 - 2*P - Q^2 + 1)/P^2*x^4 
        + (A^2*P^2*Q + 4*A*P^2 - 3*A*P*Q^2 + A*P + 4*P^2*Q - 4*P*Q 
        + 2*Q^3 - 2*Q)/P^2*x^3 + (4*A^2*P^3 - 12*A*P^2*Q + A*P*Q 
        + 16*P^3 - 8*P^2 + 8*P*Q^2 + 4*P - Q^2 + 1)/P^2*x^2 
        + (8*A*P - 12*Q)/P*x + 4/P;

    return f;

end intrinsic;

intrinsic HypG2Crv(APQ::SeqEnum) -> CrvHyp
{Returns the APQ curve}
    // return MinimalWeierstrassModel(HyperellipticCurve(Genus2Curve(APQ)));
    return MinimalWeierstrassModel(HyperellipticCurve(MinRedBinaryForm(Genus2Curve(APQ))));
end intrinsic;

intrinsic RICurve(f::RngUPolElt) -> RngUPolElt
{Spots the (isomorphic) Richelot isogenous curve}
    A<X> := Parent(f);
    C := HyperellipticCurve(f);
    K<x,y> := FunctionField(C);
    m := map<C -> Ambient(C) | [2/x, 4*y/x^3,1]>;
    C2 := Image(m);
    g := DefiningPolynomial(C2); cc := MonomialCoefficient(g, C.2^2);
    return Evaluate(DefiningPolynomial(C2), [X,0,1])/(-cc);
end intrinsic;

function KummerAndRIKummer(APQ)
    f := Genus2Curve(APQ);
    g := RICurve(Genus2Curve(APQ));
    f0,f1,f2,f3,f4,f5,f6 := Explode(Coefficients(g));
    C1 := HyperellipticCurve(f);
    C2 := HyperellipticCurve(g);
    K1 := KummerSurfaceScheme(C1);
    K2 := KummerSurfaceScheme(C2);

    K1_k := KummerSurface(Jacobian(C1));
    K2_k := KummerSurface(Jacobian(C2));

    return K1, K2, K1_k, K2_k;
end function;

intrinsic RichelotOnKummer(APQ::SeqEnum) -> MapSch
{The Richelot isogeny from the kummer surface}
    A,P,Q := Explode(APQ);
    K1, K2 := KummerAndRIKummer(APQ);
    PP<x1,x2,x3,x4> := CoordinateRing(K1);
    
    d := (-16*A^4*P^6 + A^4*P^4*Q^2 + 96*A^3*P^5*Q - 6*A^3*P^3*Q^3 + 2*A^3*P^3*Q - 
        128*A^2*P^6 - 200*A^2*P^4*Q^2 + 80*A^2*P^4 + 13*A^2*P^2*Q^4 - 14*A^2*P^2*Q^2
        + A^2*P^2 + 384*A*P^5*Q + 168*A*P^3*Q^3 - 168*A*P^3*Q - 12*A*P*Q^5 + 
        24*A*P*Q^3 - 12*A*P*Q - 256*P^6 - 240*P^4*Q^2 - 192*P^4 - 48*P^2*Q^4 + 
        96*P^2*Q^2 - 48*P^2 + 4*Q^6 - 12*Q^4 + 12*Q^2 - 4)/P^4;

    li := [
        (A*P*Q - 8*P^2 + 4*P - Q^2 + 1)/P^2*x1^2 + 8*x1*x3 + x1*x4 - 1/P*x2^2 + (-A*P + 
            2*Q)/P*x2*x3 + (P - 1)/P*x3^2,
        4*Q/P*x1^2 + (2*A*P*Q - 2*Q^2 + 2)/P^2*x1*x2 + (A^2*P^2*Q - 4*A*P^2 - 3*A*P*Q^2 
            + A*P + 4*P^2*Q + 4*P*Q + 2*Q^3 - 2*Q)/P^2*x1*x3 + (2*A*P - 2*Q)/P*x2^2 + 
            (2*A^2*P^2 - 6*A*P*Q + 8*P^2 + 4*Q^2)/P*x2*x3 + x2*x4 + (A*P - 2*Q)/P*x3^2,
        (-16*P + 4)/P*x1^2 - 4*Q/P*x1*x2 + 8/P*x1*x3 - 4*x2^2 + (A^2*P^2 - 3*A*P*Q + 
            4*P^2 + 4*P + 2*Q^2 - 2)/P*x3^2 + x3*x4,
        d*((-4*A^2*P^5 + A^2*P^4 + 1/4*A^2*P^3*Q^2 + 12*A*P^4*Q - 3*A*P^3*Q - 3/4*A*P^2*Q^3
            + 3/4*A*P^2*Q - 16*P^5 + 4*P^4 - 7*P^3*Q^2 + 2*P^3 + 2*P^2*Q^2 - 2*P^2 + 
            1/2*P*Q^4 - P*Q^2 + 1/2*P)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 
            3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 -
            13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 
            21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 
            16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 
            3/4*Q^4 - 3/4*Q^2 + 1/4)*x1^2 + (-4*A*P^4 + A*P^3 + 6*P^3*Q - 
            P^2*Q)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 
            1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4
            + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 
            21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 +
            12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 
            1/4)*x1*x2 + (A*P^3*Q + 16*P^4 - 10*P^3 - 2*P^2*Q^2 + 4*P^2)/(A^4*P^6 - 
            1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6
            + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 
            1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 
            3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 
            6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 1/4)*x1*x3 + (P^4 - 
            1/2*P^3)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 
            1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4
            + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 
            21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 +
            12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 
            1/4)*x1*x4 + (-1/4*A*P^3*Q - 2*P^4 + 1/2*P^2*Q^2 - 1/2*P^2)/(A^4*P^6 - 
            1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6
            + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 
            1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 
            3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 
            6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 1/4)*x2^2 + (A*P^4 - 
            1/2*A*P^3 - P^3*Q)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 
            3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 -
            13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 
            21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 
            16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 
            3/4*Q^4 - 3/4*Q^2 + 1/4)*x2*x3 + (1/16*A^3*P^4*Q + 1/2*A^2*P^5 - 
            5/16*A^2*P^3*Q^2 + 1/16*A^2*P^3 - 5/4*A*P^4*Q + 1/4*A*P^3*Q + 1/2*A*P^2*Q^3 
            - 1/2*A*P^2*Q + 2*P^5 - 2*P^4 + 1/2*P^3*Q^2 + 1/2*P^3 - 1/4*P^2*Q^2 + 
            1/4*P^2 - 1/4*P*Q^4 + 1/2*P*Q^2 - 1/4*P)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 
            6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2
            - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 
            24*A*P^5*Q - 21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 
            3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 
            1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 1/4)*x3^2 + (-1/2*P^4 + 1/4*P^3)/(A^4*P^6 - 
            1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6
            + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 
            1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 
            3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 
            6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 1/4)*x3*x4 - 
            1/16*P^4/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 
            1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4
            + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 
            21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 +
            12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 
            1/4)*x4^2)];

    return map<K1 -> K2 | li>;
end intrinsic;

intrinsic RichelotOnKummer(K1::Sch, K2::Sch, APQ::SeqEnum) -> MapSch
{The Richelot isogeny from the kummer surface}
    A,P,Q := Explode(APQ);
    PP<x1,x2,x3,x4> := CoordinateRing(K1);
    
    d := (-16*A^4*P^6 + A^4*P^4*Q^2 + 96*A^3*P^5*Q - 6*A^3*P^3*Q^3 + 2*A^3*P^3*Q - 
        128*A^2*P^6 - 200*A^2*P^4*Q^2 + 80*A^2*P^4 + 13*A^2*P^2*Q^4 - 14*A^2*P^2*Q^2
        + A^2*P^2 + 384*A*P^5*Q + 168*A*P^3*Q^3 - 168*A*P^3*Q - 12*A*P*Q^5 + 
        24*A*P*Q^3 - 12*A*P*Q - 256*P^6 - 240*P^4*Q^2 - 192*P^4 - 48*P^2*Q^4 + 
        96*P^2*Q^2 - 48*P^2 + 4*Q^6 - 12*Q^4 + 12*Q^2 - 4)/P^4;

    li := [
        (A*P*Q - 8*P^2 + 4*P - Q^2 + 1)/P^2*x1^2 + 8*x1*x3 + x1*x4 - 1/P*x2^2 + (-A*P + 
            2*Q)/P*x2*x3 + (P - 1)/P*x3^2,
        4*Q/P*x1^2 + (2*A*P*Q - 2*Q^2 + 2)/P^2*x1*x2 + (A^2*P^2*Q - 4*A*P^2 - 3*A*P*Q^2 
            + A*P + 4*P^2*Q + 4*P*Q + 2*Q^3 - 2*Q)/P^2*x1*x3 + (2*A*P - 2*Q)/P*x2^2 + 
            (2*A^2*P^2 - 6*A*P*Q + 8*P^2 + 4*Q^2)/P*x2*x3 + x2*x4 + (A*P - 2*Q)/P*x3^2,
        (-16*P + 4)/P*x1^2 - 4*Q/P*x1*x2 + 8/P*x1*x3 - 4*x2^2 + (A^2*P^2 - 3*A*P*Q + 
            4*P^2 + 4*P + 2*Q^2 - 2)/P*x3^2 + x3*x4,
        d*((-4*A^2*P^5 + A^2*P^4 + 1/4*A^2*P^3*Q^2 + 12*A*P^4*Q - 3*A*P^3*Q - 3/4*A*P^2*Q^3
            + 3/4*A*P^2*Q - 16*P^5 + 4*P^4 - 7*P^3*Q^2 + 2*P^3 + 2*P^2*Q^2 - 2*P^2 + 
            1/2*P*Q^4 - P*Q^2 + 1/2*P)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 
            3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 -
            13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 
            21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 
            16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 
            3/4*Q^4 - 3/4*Q^2 + 1/4)*x1^2 + (-4*A*P^4 + A*P^3 + 6*P^3*Q - 
            P^2*Q)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 
            1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4
            + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 
            21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 +
            12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 
            1/4)*x1*x2 + (A*P^3*Q + 16*P^4 - 10*P^3 - 2*P^2*Q^2 + 4*P^2)/(A^4*P^6 - 
            1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6
            + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 
            1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 
            3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 
            6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 1/4)*x1*x3 + (P^4 - 
            1/2*P^3)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 
            1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4
            + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 
            21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 +
            12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 
            1/4)*x1*x4 + (-1/4*A*P^3*Q - 2*P^4 + 1/2*P^2*Q^2 - 1/2*P^2)/(A^4*P^6 - 
            1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6
            + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 
            1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 
            3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 
            6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 1/4)*x2^2 + (A*P^4 - 
            1/2*A*P^3 - P^3*Q)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 
            3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 -
            13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 
            21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 
            16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 
            3/4*Q^4 - 3/4*Q^2 + 1/4)*x2*x3 + (1/16*A^3*P^4*Q + 1/2*A^2*P^5 - 
            5/16*A^2*P^3*Q^2 + 1/16*A^2*P^3 - 5/4*A*P^4*Q + 1/4*A*P^3*Q + 1/2*A*P^2*Q^3 
            - 1/2*A*P^2*Q + 2*P^5 - 2*P^4 + 1/2*P^3*Q^2 + 1/2*P^3 - 1/4*P^2*Q^2 + 
            1/4*P^2 - 1/4*P*Q^4 + 1/2*P*Q^2 - 1/4*P)/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 
            6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2
            - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 
            24*A*P^5*Q - 21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 
            3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 
            1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 1/4)*x3^2 + (-1/2*P^4 + 1/4*P^3)/(A^4*P^6 - 
            1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 1/8*A^3*P^3*Q + 8*A^2*P^6
            + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4 + 7/8*A^2*P^2*Q^2 - 
            1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 21/2*A*P^3*Q + 3/4*A*P*Q^5 - 
            3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 + 12*P^4 + 3*P^2*Q^4 - 
            6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 1/4)*x3*x4 - 
            1/16*P^4/(A^4*P^6 - 1/16*A^4*P^4*Q^2 - 6*A^3*P^5*Q + 3/8*A^3*P^3*Q^3 - 
            1/8*A^3*P^3*Q + 8*A^2*P^6 + 25/2*A^2*P^4*Q^2 - 5*A^2*P^4 - 13/16*A^2*P^2*Q^4
            + 7/8*A^2*P^2*Q^2 - 1/16*A^2*P^2 - 24*A*P^5*Q - 21/2*A*P^3*Q^3 + 
            21/2*A*P^3*Q + 3/4*A*P*Q^5 - 3/2*A*P*Q^3 + 3/4*A*P*Q + 16*P^6 + 15*P^4*Q^2 +
            12*P^4 + 3*P^2*Q^4 - 6*P^2*Q^2 + 3*P^2 - 1/4*Q^6 + 3/4*Q^4 - 3/4*Q^2 + 
            1/4)*x4^2)];

    return map<K1 -> K2 | li>;
end intrinsic;


intrinsic KummerIso(APQ::SeqEnum) -> MapSch
{An isomorphism iota from the Richelot isogenous kummer}
    K1, K2 := KummerAndRIKummer(APQ);

    xpx := K2.2/K2.1; xx := K2.3/K2.1; b := K2.4/K2.1;
    xpx2 := 2*xpx/xx; xx2 := 4/xx; b2 := b/xx;

    return map<K2 -> K1 | [1,xpx2,xx2, b2]>;
end intrinsic;

intrinsic KummerIso(K1::Sch, K2::Sch, APQ::SeqEnum) -> MapSch
{An isomorphism iota from the Richelot isogenous kummer}

    xpx := K2.2/K2.1; xx := K2.3/K2.1; b := K2.4/K2.1;
    xpx2 := 2*xpx/xx; xx2 := 4/xx; b2 := b/xx;

    return map<K2 -> K1 | [1,xpx2,xx2, b2]>;
end intrinsic;

intrinsic Sqrt2Kummer(APQ::SeqEnum) -> RngUPolElt, Sch, SeqEnum
{Explicit sqrt 2 mulitplication on the Kummer surface}
    K1, K2, K1_k, K2_k := KummerAndRIKummer(APQ);

    rho := RichelotOnKummer(K1,K2,APQ);
    iota := KummerIso(K1,K2,APQ);

    sqrt2 := DefiningPolynomials(rho*iota);
    return Genus2Curve(APQ), K1, sqrt2;
end intrinsic;


intrinsic Sqrt2Kummer(APQ::SeqEnum, P::SrfKumPt) -> SrfKumPt
{Explicit sqrt 2 mulitplication on the Kummer surface}
    _, K1, sqrt2 := Sqrt2Kummer(APQ);
    sqrt2 := map<K1 -> K1 | sqrt2>;
    p := Eltseq(P);
    sqrt2p := sqrt2(K1!p);
    return Parent(P)!Coordinates(sqrt2p);
end intrinsic;


intrinsic Times3Kummer(K::SrfKum) -> SeqEnum
{Multiplication by 3 on the Kummer}
    P3 := ProjectiveSpace(BaseRing(K), 3);

    two_times := [Evaluate(p, [P3.i : i in [1..4]]) : p in K`Delta];

    BB := K`BBMatrix;
    L12 := two_times cat [P3.i : i in [1..4]];
    L3 := [P3.i : i in [1..4]];
    c1 := P3.1; c2 := Evaluate(BB[1,1], L12);
    three_times := [ c1*Evaluate(BB[1,j], L12) - L3[j]*c2 : j in [1..4] ];
    three_times[1] := c1*c2;

    return three_times;
end intrinsic;

intrinsic Times3Kummer(APQ::SeqEnum) -> SeqEnum
{Multiplication by 3 on the Kummer}
    K := KummerSurface(Jacobian(HyperellipticCurve(Genus2Curve(APQ))));
    return Times3Kummer(K);
end intrinsic;


intrinsic Times2Kummer(K::SrfKum) -> SeqEnum
{Multiplication by 2 on the Kummer}
    P3 := ProjectiveSpace(BaseRing(K), 3);

    two_times := [Evaluate(p, [P3.i : i in [1..4]]) : p in K`Delta];

    return two_times;
end intrinsic;

intrinsic Times2Kummer(APQ::SeqEnum) -> SeqEnum
{Multiplication by 2 on the Kummer}
    K := KummerSurface(Jacobian(HyperellipticCurve(Genus2Curve(APQ))));
    return Times2Kummer(K);
end intrinsic;


intrinsic CoordsMulta(APQ::SeqEnum) -> RngMPolElt
{The \xi_1 coordinate of the 3 + \sqrt2 or 3 - \sqrt2  multiplication}
    K := KummerSurface(Jacobian(HyperellipticCurve(Genus2Curve(APQ))));
    P3 := ProjectiveSpace(BaseRing(K), 3);

    three := [Evaluate(p, [P3.i : i in [1..4]]) : p in Times3Kummer(APQ)];
    _,_,sqrt2 := Sqrt2Kummer(APQ);
    sqrt2 := [Evaluate(p, [P3.i : i in [1..4]]) : p in sqrt2];

    BB := K`BBMatrix;

    L12 := three cat sqrt2;
    return Evaluate(BB[1,1], L12), Evaluate(BB[2,2], L12), Evaluate(BB[3,3], L12);
    
end intrinsic;



function UserChoicePlys(ff)
    fs := [f : f in ff | Degree(f[1]) le 24];

    read choice, "We need your help, choose plys from " cat Sprint(fs);
    choice := eval choice;

    g := &*[fs[i][1] : i in choice[1]];
    h := &*[fs[i][1] : i in choice[2]];

    return g, h;
end function;

function Deg24Factors(ff)
    ff2 := &cat[ [f[1] : i in [1..f[2]]] : f in ff | Degree(f[1]) le 24];

    to24 := [ss : ss in Subsets({i : i in [1..#ff]}) | &+([0] cat [Degree(ff2[i]) : i in ss]) eq 24];

    assert #to24 eq 2;
    g := &*[ff[i] : i in to24[1]];
    h := &*[ff[i] : i in to24[2]];

    return g, h;
end function;

intrinsic Poly24Alg(APQ::SeqEnum) -> SeqEnum
{Computes four polynomials, the roots of the ...}
    f := Genus2Curve(APQ);  
    C := HyperellipticCurve(f); KK := BaseRing(C);
    J := Jacobian(C); K := KummerSurface(J);

    three := Times3Kummer(APQ); three := [three[2]/three[1], three[3]/three[1], three[4]/three[1]];
    _,_,sqrt2 := Sqrt2Kummer(APQ);  sqrt2 := [sqrt2[2]/sqrt2[1], sqrt2[3]/sqrt2[1], sqrt2[4]/sqrt2[1]];

    P3 := PolynomialRing(KK, 4);
    fs := [P3!Numerator(three[i] - sqrt2[i]) : i in [1..3]];
    fK := P3!DefiningPolynomial(K);

    A3<x,y,z> := PolynomialRing(KK, 3);  _<X> := PolynomialRing(KK);


    fs := [Evaluate(ff, [1,x,y,z]) : ff in fs];
    fK := Evaluate(fK, [1,x,y,z]);

    r1 := Resultant(fK, fs[1], 3);
    r2 := Resultant(fK, fs[2], 3);

    //the sum
    vprint TalkToMe, 1 : "First resultant";
    plys := Resultant(r1, r2, 2);
    ff := Factorisation(plys);


    gh1 := [g[1] : g in ff | Degree(g[1]) eq 24];
    if not #gh1 in [1,2] then
        vprint TalkToMe, 1 : "Yeah, sorry, both Galois representations are reducible. Need to fix this.";
        return [[],[]];
    end if;

    gh1 := [Evaluate(g, [X,0,0]) : g in gh1];

    //the product
    vprint TalkToMe, 1 : "Second resultant";
    plys := Resultant(r1, r2, 1);
    ff := Factorisation(plys);

    gh2 := [g[1] : g in ff | Degree(g[1]) eq 24];
    if not #gh2 eq #gh1 then
        "Should never get here.";
        return [[],[]];
    end if;

    gh2 := [Evaluate(g, [0,X,0]) : g in gh2];

    // Now want to order the polynomials correctly, so that they define the same field
    p := NextPrime(100);
    good := false;
    while not good do
        p := NextPrime(p);
        if not p in BadPrimes(IntegralModel(C)) then
            good := true;
        end if;
    end while;
    Fp := GF(p);
    P := PolynomialRing(Fp);
    Fq := SplittingField(P!gh1[1]);
    rts := [r[1] : r in Roots((P!gh1[1]), Fq)];
    rts2 := [r[1] : r in Roots((P!gh2[1]), Fq)];
    Kum_p := BaseChange(K, Fq);
    pts_set := &join{&join{{pt : pt in Points(Kum_p, [1,a,b]) | 7*pt eq Kum_p!0} : b in rts2} : a in rts};

    // Should have 24 points = (49-1)/2, so swap them if not
    if #pts_set eq 24 then     
        return [gh1,gh2];
    else
        return [gh1,[gh2[2],gh2[1]]];
    end if;
    
end intrinsic;

