AttachSpec("../spec");

//////////////////////////////////////////////////
// Proof of Lemma 3.6.5
//////////////////////////////////////////////////

for N in [67..2500] do
    k := Valuation(N,2);
    M := N div 2^k;

    mu := M^2*(&*([1] cat [1 + 1/p : p in PrimeDivisors(M)]));
    
    pg_bound := 1/2*(1/480*N*(N-1)*(N-23) - 1/4*2^(2*k-1)*1/2*mu - 1);

    assert pg_bound ge 3;
end for;

// the values of p_g for N \leq 24 are recorded in the tables (see the
// tables directory for how to compute these) so it suffices to check
// claim for N > 25
for N in [25, 67] do
    for r in [r : r in [1..N] | GCD(N,r) eq 1] do
        assert GeometricGenusWNr(N, r) ge 3;
    end for;
end for;


//////////////////////////////////////////////////
// Proof of Lemma 3.6.6
//////////////////////////////////////////////////
for N in [106..6000] do
    k := Valuation(N,2);
    M := N div 2^k;

    mu := M^2*(&*([1] cat [1 + 1/p : p in PrimeDivisors(M)]));
    
    K2_bound := 1/2*(1/60*N*(N-1)*(N-30) - 3*2^(2*k-1)*1/2*mu - 6);

    assert K2_bound gt 0;
end for;
