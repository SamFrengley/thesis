/*
    Using the MW Sieve to prove the last part of Theorem 4.5.1.
*/

load "../functions/mwsieve.m";
E30 := EllipticCurve([1,-1,0,-273176601587417,-1741818799948905109620]);
E48 := EllipticCurve([0,0,1, 468240736152891010, -148374586624464876247316957]);
j30 := jInvariant(E30);
j48 := jInvariant(E48);

// ----------------
// X+(ns3, ns5)
// ----------------

E := eval Read("../models/eqns/X_ns3+-ns5+.m"); assert CremonaReference(E) eq "225a1";
j := eval Read("../models/jmaps/jX_ns3+-ns5+.m");

_<u,v> := FunctionField(E);
d := 3*(3*u^2 - 2*u - 4*v + 3)*(9*u^2 - 10*u + 5);

gens, pts, bound := MWSieveInfo(E, d : B:=10000);

for pt in pts do
    jp := j(pt);
    if jp[2] ne 0 then 
        jp := jp[1]/jp[2];
        assert HasComplexMultiplication(EllipticCurveWithjInvariant(jp)) or jp eq j30;
    end if;
end for;

can_height := bound[1]^2*CanonicalHeight(gens[1]);
assert can_height gt 10^(15);

// ----------------
// X+(ns3, s5)
// ----------------

E := eval Read("../models/eqns/X_ns3+-s5+.m"); assert CremonaReference(E) eq "225a1";
j := eval Read("../models/jmaps/jX_ns3+-s5+.m");

_<u,v> := FunctionField(E);
d := -(3*v + 3*u^3 - 12*u^2 + 12*u - 6);

gens, pts, bound := MWSieveInfo(E, d : B:=10000);

for pt in pts do
    jp := j(pt);
    if jp[2] ne 0 then 
        jp := jp[1]/jp[2];
        assert HasComplexMultiplication(EllipticCurveWithjInvariant(jp));
    end if;
end for;

can_height := bound[1]^2*CanonicalHeight(gens[1]);
assert can_height gt 10^(15);

// ----------------
// X+(ns3, ns7)
// ----------------

E := eval Read("../models/eqns/X_ns3+-ns7+.m"); assert CremonaReference(E) eq "441b1";
j := eval Read("../models/jmaps/jX_ns3+-ns7+.m");

_<u,v> := FunctionField(E);
d := 4*u^3 + 21*u^2 - 42*u + 49;

gens, pts, bound := MWSieveInfo(E, d : B:=5000);

for pt in pts do
    jp := j(pt);
    if jp[2] ne 0 then 
        jp := jp[1]/jp[2];
        assert HasComplexMultiplication(EllipticCurveWithjInvariant(jp));
    end if;
end for;

can_height := bound[2]^2*CanonicalHeight(gens[2]);
assert can_height gt 10^(15);

// ----------------
// X+(ns2, ns11)
// ----------------

E := eval Read("../models/eqns/X_ns11+.m"); assert CremonaReference(E) eq "121b1";
j := eval Read("../models/jmaps/jX_ns11+.m");

_<u,v> := FunctionField(E);
d := (-3*u^2 + 24*u - 37)*v - u^4 + 5*u^3 + 18*u^2 - 95*u + 94;

gens, pts, bound := MWSieveInfo(E, d : B:=10000);

for pt in pts do
    jp := j(pt);
    if jp[2] ne 0 then 
        jp := jp[1]/jp[2];
        assert HasComplexMultiplication(EllipticCurveWithjInvariant(jp));
    end if;
end for;

can_height := bound[1]^2*CanonicalHeight(gens[1]);
assert can_height gt 10^(15);

// ----------------
// X(ns8, ns3)/delta
// ----------------

E := eval Read("../models/eqns/X_ns8+-ns3+.m"); assert CremonaReference(E) eq "576a1";
j := eval Read("../models/jmaps/jX_ns8+-ns3+.m");

_<u,v> := FunctionField(E);
d := 3*(3*u^3 + 3*u^2 + 2*u*v - 8*u + 8);

gens, pts, bound := MWSieveInfo(E, d : B:=10000);

for pt in pts do
    jp := j(pt);
    if jp[2] ne 0 then 
        jp := jp[1]/jp[2];
        assert HasComplexMultiplication(EllipticCurveWithjInvariant(jp)) or jp eq j48;
    end if;
end for;

can_height := bound[2]^2*CanonicalHeight(gens[2]);
assert can_height gt 10^(15);
