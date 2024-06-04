/*
    Suppose that there exists a morphism phi : C -> E defined over Q. Let P1,..., Pl be generators for 
    E(Q). This function returns a set S of points together with a bound N = (N1,...,Nl) 
    such that if phi^(-1)(P) \in C(Q) then either P \in S or we have P = a1*P1 + ... + al*Pl where one of 
    the ai is outside the range (-Ni, Ni).

    Notable functions are:

    MWSieve(E, d : B:=5000, vSmooth:=7)
        // Given an elliptic curve E/Q given by Weierstrass equation f(u,v) = 0 
        // consider the double cover C/Q : {f(u,v) = 0 and w^2 = d(u,v)} where d 
        // is a rational function on E. Performs the MW sieve on E for primes <=B.
        // Returns: [P1,...,Pn], generators of E(Q), allowable tuples (a1,..,an) such 
        // that a1*P1 + ... + an*Pn might have a preimage on C(Q). These are ordered by 
        // max(|a_i|).
        // vSmooth is a bound on the smoothness of the order of P_i modulo p - if there
        // is a prime factor bigger than vSmooth we ignore the prime. In practice this seems to 
        // speed up the computation.

    MWSieveInfo(E, d : B:=5000, vSmooth:=7, Verbose:=false)
        // Performs the MWSieve but outputs more info. The inputs are the same as MWSieve.
        // The returns are the generators of E/Q, the points of small height which might lift,
        // and the next smallest point (not in the set of small height points found).


    No efficiency claims are being made. This is rather naive.
*/

function SortFunction(a, b)
    // Orders tuples by the maximum absolute value of their entries.
    return Max([Abs(ai) : ai in a]) - Max([Abs(bi) : bi in b]);
end function;

function CorrectSignModN(a, N)
    // Given a = (a_1,...,a_n) \in \prod_{i=1}^n Z/N_iZ chooses a lift
    // a' \in Z^n such that each a'_i \in (-N/2, N/2].

    ret := [];
    for i in [1..#N] do
        ai := a[i];
        if ai gt N[i]/2 then
            ai := ai - N[i];
        end if;
        Append(~ret, ai);
    end for;

    return ret;
end function;

function CreateProductSet(N_p)
    // Creates the set \prod_{N \in N_p} Z/N  

    number := &*[N_pi : N_pi in N_p];
    up_to_prods := [1] cat [&*N_p[1..i] : i in [1..#N_p-1]];
    S := { };

    for n in [1..number] do
        l := [(Floor(n/up_to_prods[i]) + 1) mod N_p[i] : i in [1..#N_p]];
        Include(~S, l);
    end for;

    return S;
end function;

function CanLift(P_p, N_p, d_p) 
    // Takes entries P_p = (P_1,...,P_n) the generators of E reduced modulo p, 
    // N_p = (Ord(P_i))_{i=1}^n their orders, and d_p the function d reduced modulo p.
    // Returns the set of a_p \in \prod_{N \in N_p} Z/N such that there is a point
    // on \tilde{C} : { f(u,v) = 0 and w^2 - d(u,v) = 0} mod p 
    // above Q = \sum_{a_i \in a_p} a_i P_i. 

    curlyP_p := CreateProductSet(N_p);
    A_p := { [0 : i in N_p] };

    for a_p in [a : a in curlyP_p] do
        d_here := Evaluate(d_p, &+[a_p[i]*P_p[i] : i in [1..#N_p]]);

        if Type(d_here) ne FldFinElt then //this is at the poles of d
            Include(~A_p, a_p);
        elif IsSquare(d_here) then
            Include(~A_p, a_p);
        end if;
    end for;

    return A_p;
end function;

function PatchTogether(A, N, A_p, N_p)
    // Takes entries N = (N_1,...,N_n) and N_p = (N_{p,1},...,N_{p,n}) and
    // A, A_p are allowable congruence classes modulo N and N_p respectively. 
    // Returns the allowable congruence calsses modulo M = (LCM(N_i, N_{p,i}))
    // and M.
    M := [LCM(N[i], N_p[i]) : i in [1..#N]];

    S := {};

    for a_p in A_p do
        for a in A do
            b := [];
            flag := true; i := 1;
            while flag do
                b_i := CRT([a_p[i], a[i]], [N_p[i], N[i]]);
                i := i + 1;
                flag := b_i ne -1; ;
                if flag then Append(~b, b_i); end if;
                flag := flag and (i le #N);
            end while;

            if #b eq #N then
                Include(~S, b);
            end if;
        end for;
    end for;

    return S, M;
end function;

function MWSieve(E, d : B:=5000, vSmooth:=7)
    // Given an elliptic curve E/Q given by Weierstrass equation f(u,v) = 0 
    // consider the double cover C/Q : {f(u,v) = 0 and w^2 = d(u,v)} where d 
    // is a rational function on E. Performs the MW sieve on E for primes <=B.
    // Returns: [P1,...,Pn], generators of E(Q), allowable tuples (a1,..,an) such 
    // that a1*P1 + ... + an*Pn might have a preimage on C(Q). These are ordered by 
    // max(|a_i|).
    // vSmooth is a bound on the smoothness of the order of P_i modulo p - if there
    // is a prime factor bigger than vSmooth we ignore the prime. In practice this seems to 
    // speed up the computation.

    A1 := AffineSpace(Rationals(), 1);
    dd := map<E -> A1 | [d]>; dd := DefiningPolynomials(dd)[1];

    small_primes := [p : p in PrimesInInterval(1, B) | not p in BadPrimes(E)];

    P := Generators(E); P := [Coordinates(Pi) : Pi in P];
    P := [[Integers()!(cc*LCM([Denominator(c) : c in Pi])) : cc in Pi] : Pi in P];

    aE := aInvariants(E);
    aE := [Integers()!aa : aa in aE];

    A := { [0 : i in P] };
    N := [1 : i in P];

    for p in small_primes do
        F_p := GF(p);
        E_p := EllipticCurve([F_p!aa : aa in aE]); //reduce E mod p
        P_p := [E_p![F_p!c : c in Pi] : Pi in P];
        N_p := [Order(Pi) : Pi in P_p];
        if Max([q[1] : q in Factorisation(&*N_p)]) le vSmooth then 
            _<x,y> := FunctionField(E_p);
            d_p := Evaluate(dd, [x,y,1]);

            A_p := CanLift(P_p, N_p, d_p);

            A, N := PatchTogether(A, N, A_p, N_p);
        end if;
    end for;

    A := [CorrectSignModN(a, N) : a in A];

    return [E!Q : Q in P], Sort([a : a in A], SortFunction);
end function;

function MWSieveInfo(E, d : B:=5000, vSmooth:=7, Verbose:=false)
    // Performs the MWSieve but outputs more info. The inputs are the same as MWSieve.
    // The returns are the generators of E/Q, the points of small height which might lift,
    // and the next smallest point (not in the set of small height points found).

    P, A := MWSieve(E, d : B:=B, vSmooth:=vSmooth);

    flag := true;
    small_height := [];
    i := 1;
    while flag do
        Append(~small_height, &+[A[i][j]*P[j] : j in [1..#P]]);
        i := i + 1;
        if &+[Abs(a) : a in A[i]] gt 10^3 then
            flag := false;
        end if;
    end while;

    if Verbose then
        Sprintf("Unless there is a point further away than %o (generators are %o) then the only points which may have a nonempty fibre are: \n%o", A[i], P, small_height);
    end if;

    return P, small_height, A[i];    
end function;
