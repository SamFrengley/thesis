/*
    Functions useful for manipulating subgroups of GL2(Z/N) for the 
    purpose of determining QTT congruences.

    Notable functions are

    GL2GivesQTTCongruence(Hplus)
        // Given a subgroup H+ \subset GL2(Z/N) returns two values. The first is a 
        // boolean, true if and only if there is a subgroup and g satisfying (1) in
        // Lemma 3.1. The second is a list of pairs <H,gg> where gg is the set of all 
        // g for which (1) holds for the pair H+,H. 

    IsInducedCongruence(H1plus, H1, H2plus, H2)
        // Given two pairs of subgroups (H1,H1+) and (H2,H2+) of GL2(Z/N) this function 
        // determines if there exists a g \in GL2(Z/N) such that conjugating H2 and H2+ by 
        // g we have: 
        // g^(-1) H2+ g \subset H1+ 
        // g^(-1) H2 g = H1 \cap (g^(-1) H2+ g)
        // the function returns a boolean, and g.

    GL2ModDeltaGroup(H1, H1p, H2, H2p)
        // Given pairs of subgroups H1,H1+ \subset GL2(Z/N1) and H2,H2+ \subset GL2(Z/N2) satisfying 
        // Lemma 3.1 equation (1) such that (N1, N2) = 1 returns the subgroup corresponding to the 
        // modular curve X(H1, H2)/delta in GL2(Z/N1*N2). 
*/

// Attach Rouse--Sutherland--Zureick-Brown's (super useful!!) package.
Attach("../RSZ21/gl2.m"); 

function RZBLoad()
    // Loads the 1208 subgroups of GL2(Z_2) given by Rouse--Zureick-Brown which occur as
    // the 2-adic image of the Galois rep attached to a non-CM elliptic curve /Q.
    // Altered from a function of the same name in RSZ21.

    rzbdata := [Split(s,":") : s in Split(Read("../RSZ21/rzbdata.txt"))];
    rzbdata := [<r[1],GL2FromGenerators(StringToInteger(r[2]),eval(r[3]))> : r in rzbdata];
    rzbnonCM := eval Read("../RSZ21/rzbimages.txt");
    return [grp : grp in rzbdata | grp[1] in rzbnonCM];
end function;

function ComputeTheg(Hplus, H)
    // Given H,H+ \susbet GL2(Z/N) this computes the g for which equation (1)
    // in Lemma 3.1 holds.
    N := #(BaseRing(Hplus));
    G := GL(2, Integers(N));
    assert exists(h){h : h in Hplus | not h in H};

    glist := [];
    for g in [g : g in Centralizer(G, H) | 2*Trace(g) eq 0] do
        if g^(-1)*h*g eq -h then
            truth := true;
        else 
            truth := false;
        end if;

        if truth then
            Append(~glist, G!g);
        end if;
    end for;

    return glist;
end function;

function GL2GivesQTTCongruence(Hplus)
    // Given a subgroup H+ \subset GL2(Z/N) returns two values. The first is a 
    // boolean, true if and only if there is a subgroup and g satisfying (1) in
    // Lemma 3.1. The second is a list of pairs <H,gg> where gg is the set of all 
    // g for which (1) holds for the pair H+,H. 

    N := #(BaseRing(Hplus));
    G := GL(2, Integers(N));

    Hprime := sub<Hplus | [h : h in Hplus | 2*Trace(h) ne 0]>;

    if IsOdd(Index(Hplus, Hprime)) then
        return false, [];
    end if;

    if Index(Hplus, Hprime) gt 2 then
        candidates := [GL2MinimizeGenerators(H) : H in IntermediateSubgroups(Hplus, Hprime) | Index(Hplus, H) eq 2];
    else 
        candidates := [GL2MinimizeGenerators(Hprime)];
    end if;
    
    congs := [];

    for H in candidates do
        g_for_H := ComputeTheg(Hplus, H);

        if #g_for_H gt 0 then
            Append(~congs, <H, g_for_H>);
        end if;
    end for;

    return (#congs gt 0), congs;
end function; 

function IsInducedCongruence(H1plus, H1, H2plus, H2)
    // Given two pairs of subgroups (H1,H1+) and (H2,H2+) of GL2(Z/N) this function 
    // determines if there exists a g \in GL2(Z/N) such that conjugating H2 and H2+ by 
    // g we have: 
    // g^(-1) H2+ g \subset H1+ 
    // g^(-1) H2 g = H1 \cap (g^(-1) H2+ g)
    // the function returns a boolean, and g.

    N := #(BaseRing(H1plus));
    G := GL(2, Integers(N));

    flag := exists(H2plus_prime){H : H in Class(G, H2plus) | H subset H1plus};
    if not flag then return flag; end if;

    _, g := IsConjugate(G, H2plus, H2plus_prime);
    H2_prime := Conjugate(H2, g);

    N_G_H2plus_prime := Normaliser(G, H2plus_prime);

    flag, g2 := IsConjugate(N_G_H2plus_prime, H2_prime, (H2plus_prime meet H1));
    if not flag then return flag; end if;

    assert H2plus_prime eq Conjugate(H2plus, g2*g);
    assert (H2plus_prime meet H1) eq Conjugate(H2, g2*g);

    return flag, g2*g;
end function;

function GL2Patch(H1, H2)
    // Given subgroups H1 \subset GL2(Z/N1) and H2 \subset GL2(Z/N2) such that (N1, N2) = 1
    // returns the subgroup corresponding to the modular curve X(H1, H2) in GL2(Z/N1*N2). 
    // Equivalently this is the subgroup given by the intersection of their inverse images
    // under GL2(Z/N1*N2) -> GL2(Z/N1), resp GL2(Z/N1*N2) -> GL2(Z/N2).

    N1 := #(BaseRing(H1));
    N2 := #(BaseRing(H2));
    assert GCD(N1, N2) eq 1;
    N := N1*N2;
    return GL2Lift(H1, N) meet GL2Lift(H2, N);
end function;

function GL2ModDeltaGroup(H1, H1p, H2, H2p)
    // Given pairs of subgroups H1,H1+ \subset GL2(Z/N1) and H2,H2+ \subset GL2(Z/N2) satisfying 
    // Lemma 3.1 equation (1) such that (N1, N2) = 1 returns the subgroup corresponding to the 
    // modular curve X(H1, H2)/delta in GL2(Z/N1*N2). 

    N1 := #(BaseRing(H1));
    N2 := #(BaseRing(H2));
    assert GCD(N1, N2) eq 1;
    assert Index(H1p, H1) eq 2 and Index(H2p, H2) eq 2;
    Hplus := GL2Patch(H1p, H2p);
    H := GL2Patch(H1, H2);
    assert exists(Hdelta){h : h in IntermediateSubgroups(Hplus, H) | Index(Hplus, h) eq 2 and 
                                                                        GL2Project(h, N1) eq H1p and
                                                                        GL2Project(h, N2) eq H2p};
                                                                    
    return Hdelta;
end function;

