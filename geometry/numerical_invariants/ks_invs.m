/**************************************************
NUMERICAL INVARIANTS FROM THE KS PAPER
**************************************************/
function m(N)
    return N^3/2*&*[1 - 1/p^2 : p in PrimeDivisors(N)];
end function;

function r_0(N)
    return ClassNumber(-4*N^2);
end function;

function r_1(N)
    return ClassNumber(-3*N^2);
end function;

function r_inf(N)
    return 1/2*(&+[EulerPhi(d)*EulerPhi(N div d) : d in Divisors(N)]) - 1/2*EulerPhi(N);
end function;

function s_11e(N,r)
    if N mod 3 ne 0 then
        return (1/2)*r_1(N);
    else
        assert r mod 3 ne 0;
        if r mod 3 eq 1 then
            return r_1(N);
        else
            return 0;
        end if;
    end if;
end function;

function LeastPositiveResidue(a, b)
    l := [1..a];
    x := Min([u : u in l | b mod a eq u]);
    return x;
end function;

function bbL(n, q)
    qq := LeastPositiveResidue(n,q);
    return #HJContinuedFraction(n/qq);
end function;

function bbL_infe(N,r)
    return &+[ EulerPhi(d)/2 * &+[ bbL(N, r*n^2*d) : n in [1..(N div d)] | GCD(n, N div d) eq 1 ] : d in Divisors(N) | d ne N];
end function;

function bbL_e(N,r)
    return r_0(N) + 2*r_1(N) - s_11e(N,r) + bbL_infe(N,r);
end function;

function sawtooth(x)
    if Denominator(x) eq 1 then
        return 0;
    else
        return x - Floor(x) - 1/2;
    end if;
end function;

function bbS(q,n)
    return &+[ sawtooth(k/n)*sawtooth(k*q/n) : k in [1..n-1]];
end function;

function bbS_infe(N,r)
    return &+[ EulerPhi(d)/2 * &+[ bbS(r*n^2, N div d) : n in [1..(N div d)] | GCD(n, N div d) eq 1] : d in Divisors(N) | d ne N];
end function;

function R_infe(N,r)
    return 12*bbS_infe(N,r) + bbL_infe(N,r) + r_inf(N);
end function;

function bbS_e(N,r)
    return 1/18*(2*s_11e(N,r) - r_1(N)) + bbS_infe(N,r);
end function;

function IntPart(x)
    return x - Floor(x);
end function;

intrinsic Recursion_ai(N::RngIntElt, d::RngIntElt,q::RngIntElt) -> SeqEnum[RngIntElt]
{
The multiplicities of the components of the singularity of type d/q in 
the divisor j*(infty) resp (j')*(infty)
}
    assert exists(qq){qq : qq in [1..d-1] | Integers(d)!qq eq Integers(d)!q};
    c := HJContinuedFraction(d/qq);
    
    M := [ [1] cat [0 : i in [1..#c+1]] ];
    for i in [1..#c] do
        Append(~M, [0 : j in [1..i-1]] cat [-1,c[i],-1] cat [0 : j in [1..(#c-i)]]);
    end for;
    Append(~M, [0 : i in [1..#c+1]] cat [1]);

    M := Matrix(Integers(), #c+2, #c+2, M);
    M := Transpose(M);

    S := [N] cat [0 : i in [1..#c+1]];
    S := Matrix(Integers(),1,#c+2,S);

    x := Eltseq(Solution(M,S));

    S := [0 : i in [1..#c+1]] cat [N];
    S := Matrix(Integers(),1,#c+2,S);

    y := Eltseq(Solution(M,S));
    
    return x,y;
end intrinsic;


intrinsic varrho(N::RngIntElt,m::RngIntElt) -> RngIntElt
{
The correction term varrho(N,m) in the formular for K.Ftil_m
}
  ST := [];
  kks := [];
    for n in Divisors(m) do
        assert exists(q){q : q in {1..N-1} | Integers(N)!m eq Integers(N)!(q*n^2)};
        a,ap := Recursion_ai(N,N,q);
        ap2 := [10^(-10)] cat [a : a in ap[2..#ap]];
        
        assert exists(k){k : k in [2..#a] | (a[k-1]/ap2[k-1] gt m/n^2) and (a[k]/ap2[k] le m/n^2)};
        Append(~kks, <q,k-1>);

        M := Matrix(Integers(), 2, 2, [a[k], ap[k], a[k-1], ap[k-1]]); 
        S := Matrix(Integers(), 1, 2, [m div n, n]);
        x := Solution(M,S);
        Append(~ST, Eltseq(x));
    end for;
    
    ret := &+[EulerPhi(GCD(st[1], st[2]))*((st[1] + st[2])/GCD(st[1], st[2]) - 1) : st in ST];

    return Integers()!ret, ST, kks;              
end intrinsic;

/**************************************************
INVARIANTS OF THE SURFACE \tilde{Z}_Nr
**************************************************/

intrinsic CinfSquared(N::RngIntElt,r::RngIntElt) -> RngIntElt
{Self intersection of the curve C_(infty,1) on til(Z)_Nr, note the 1/2}
    ret :=  -&+[1/2* EulerPhi(GCD(v,N))* IntPart(v^2*r/(N*GCD(N,v))) : v in [1..N-1]];
    return Integers()!ret;
end intrinsic;


intrinsic X1Genus(N::RngIntElt) -> RngIntElt
{The genus of the modular curve X_1(N)}
    if N le 5 then
        return 0;
    else    
        ret := 1;
        ret +:= N^2/24*(&*[1-1/p^2 : p in PrimeDivisors(N)]);
        ret +:= -1/4*(&+[EulerPhi(v)*EulerPhi(N div v) : v in Divisors(N)]);
        return Integers()!ret;
    end if;
end intrinsic;


intrinsic GeometricGenusZNr(N::RngIntElt,r::RngIntElt) -> RngIntElt
{This is the code version of Kani--Schanz Theorem 2}
    ret := m(N)*(N-12)/(144*N);
    ret +:= EulerPhi(N)/8;
    ret +:= r_0(N)/8;
    ret +:= 2*r_1(N)/9;
    ret +:= -s_11e(N,r)/9;
    ret +:= r_inf(N)/3;
    ret +:= bbL_infe(N,r)/12;
    ret +:= -R_infe(N,r)/12;
    ret +:= -1;
    return Integers()!ret;
end intrinsic;


intrinsic KZtilSquared(N::RngIntElt,r::RngIntElt) -> RngIntElt
{This is code version of Kani--Schanz Theorem 2.6}
    ret := m(N)*(N-12)/(18*N);
    ret +:= EulerPhi(N);
    ret +:= -s_11e(N,r)/3;
    ret +:= 3*r_inf(N);
    ret +:= -R_infe(N,r);
    return Integers()!ret;
end intrinsic;


intrinsic Cinf_KZtil(N::RngIntElt,r::RngIntElt) -> RngIntElt
{Compute the intersection ~C_1 dot Ktil via Prop 2.5(d)}
    return 2*X1Genus(N) - 2 - CinfSquared(N,r);
end intrinsic;


