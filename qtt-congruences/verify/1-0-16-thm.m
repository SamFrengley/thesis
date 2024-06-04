/*
    Proof that the curves E_{30}, E_{32a}, E_{32b}, E_{48} do not admit M congruences 
    over Q where M > N.
*/

function TrFrobDiff(E1, E2, p)
    return TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p);
end function;

E30 := EllipticCurve([1,-1,0,-273176601587417,-1741818799948905109620]);
E32a := EllipticCurve([0,0,1,-3014658660, 150916472601529]);
E32b := EllipticCurve([0,0,1,-26668521183591560,-1676372876256754456631251]);
E48 := EllipticCurve([0,0,1, 468240736152891010, -148374586624464876247316957]);

E30_d := QuadraticTwist(E30, -214663);
E32a_d := QuadraticTwist(E32a, Discriminant(E32a));
E32b_d := QuadraticTwist(E32b, Discriminant(E32b));
E48_d := QuadraticTwist(E48, Discriminant(E48));

// The primes <1000 which are good for all the above E and don't divide N
bad := Conductor(E30)*Conductor(E32a)*Conductor(E32b)*Conductor(E48)*30*32*48;
ll := [ell : ell in PrimesInInterval(2, 1000) | (bad mod ell) ne 0];

// Computing the GCD of the differences in the trace of frobenius to conclude
assert GCD([TrFrobDiff(E30, E30_d, ell) : ell in ll]) eq 30;
assert GCD([TrFrobDiff(E32a, E32a_d, ell) : ell in ll]) eq 32;
assert GCD([TrFrobDiff(E32b, E32b_d, ell) : ell in ll]) eq 32;
assert GCD([TrFrobDiff(E48, E48_d, ell) : ell in ll]) eq 48;