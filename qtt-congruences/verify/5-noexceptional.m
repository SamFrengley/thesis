/*
    Checking that the curves claimed to have no exceptional points in 
    Section 4.5, have no exceptional points.
*/

// ----------------
// X(G11, s3+)
// ----------------

E := eval Read("../models/eqns/X_G11-s3+.m"); assert CremonaReference(E) eq "48a1";
j := eval Read("../models/jmaps/jX_G11-s3+.m");

assert Rank(E) eq 0;
T, t := TorsionSubgroup(E);
T := {t(P) : P in T};

cm_discriminants := {};

for pt in T do
    jp := j(pt); 
    if jp[2] ne 0 then 
        jp := jp[1]/jp[2];
        has_cm, D := HasComplexMultiplication(EllipticCurveWithjInvariant(jp));
        assert has_cm;
        Include(~cm_discriminants, D);
    end if;
end for;

assert cm_discriminants eq {};

// ----------------
// X(s4+, ns3+)
// ----------------

E := eval Read("../models/eqns/X_s4+-ns3+.m"); assert CremonaReference(E) eq "36a1";
j := eval Read("../models/jmaps/jX_s4+-ns3+.m");

assert Rank(E) eq 0;
T, t := TorsionSubgroup(E);
T := {t(P) : P in T};

cm_discriminants := {};

for pt in T do
    jp := j(pt); 
    if jp[2] ne 0 then 
        jp := jp[1]/jp[2];
        has_cm, D := HasComplexMultiplication(EllipticCurveWithjInvariant(jp));
        assert has_cm;
        Include(~cm_discriminants, D);
    end if;
end for;

assert cm_discriminants eq {-4, -7};
