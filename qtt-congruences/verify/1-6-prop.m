/*
    The computations proving Proposition 4.1.6.
*/

QQ := Rationals();
Attach("../RSZ21/gl2.m");
load "../functions/qttcong.m";

// Load Rouse--Zureick-Brown's 1208 images.
rzbdata := RZBLoad(); assert #rzbdata eq 1208;

// ----------------
// LEVEL 4
// ----------------

// Allowable level 4 images
level4grps := [<grp[1], GL2Lift(grp[2], 4)> : grp in rzbdata | GL2Level(grp[2]) le 4]; 
level4grps := level4grps cat [<grp[1], GL2Project(grp[2], 4)> : grp in rzbdata | GL2Level(grp[2]) gt 4]; 

// Load level 4 data from Table 4.1
QTT41_congs := eval Read("table3-1data/QTT4-1congruencedata.txt"); 
QTT43_congs := eval Read("table3-1data/QTT4-3congruencedata.txt");


for grp in level4grps do
    flag, H_and_g := GL2GivesQTTCongruence(grp[2]);
    for data in H_and_g do
        for power in {Integers()!Determinant(g) mod 8 : g in data[2]} do // if a pair (grp, data[1]) gives a QTT (4, power)-congruence
                                                                         // then consider it
            if power eq 1 then
                // ``IsInducedCongruence'' checks that we can conjuagate (i.e., choose basis for E[4]) such that grp \subset H+ and
                // grp \cap H = data[1]
                assert exists(HH){HH : HH in QTT41_congs | IsInducedCongruence(HH[2], HH[3], grp[2], data[1])};
            
            elif power eq 3 then
                assert exists(HH){HH : HH in QTT43_congs | IsInducedCongruence(HH[2], HH[3], grp[2], data[1])};
            end if;

        end for;
    end for;
end for;

// ----------------
// LEVEL 8
// ----------------

// Allowable level 4 images
level8grps := [<grp[1], GL2Lift(grp[2], 8)> : grp in rzbdata | GL2Level(grp[2]) le 8]; 
level8grps := level8grps cat [<grp[1], GL2Project(grp[2], 8)> : grp in rzbdata | GL2Level(grp[2]) gt 8]; 

// Load level 8 data from Table 4.1
QTT81_congs := eval Read("table3-1data/QTT8-1congruencedata.txt");
QTT83_congs := eval Read("table3-1data/QTT8-3congruencedata.txt");
QTT85_congs := eval Read("table3-1data/QTT8-5congruencedata.txt"); 
QTT87_congs := eval Read("table3-1data/QTT8-7congruencedata.txt");

all_QTT8_congs := {};  //this records the 2-adic images which give rise to QTT 8-congruences

for grp in level8grps do
    flag, H_and_g := GL2GivesQTTCongruence(grp[2]);
    for data in H_and_g do
        for power in {Integers()!Determinant(g) mod 8 : g in data[2]} do

            if power eq 1 then
                assert exists(HH){HH : HH in QTT81_congs | IsInducedCongruence(HH[2], HH[3], grp[2], data[1])};
            
            elif power eq 3 then
                assert exists(HH){HH : HH in QTT83_congs | IsInducedCongruence(HH[2], HH[3], grp[2], data[1])};

            elif power eq 5 then
                assert exists(HH){HH : HH in QTT85_congs | IsInducedCongruence(HH[2], HH[3], grp[2], data[1])};

            elif power eq 7 then
                assert exists(HH){HH : HH in QTT87_congs | IsInducedCongruence(HH[2], HH[3], grp[2], data[1])};

            end if;

            Include(~all_QTT8_congs, grp[1]); 
        end for;
    end for;
end for;


// ----------------
// LEVEL 16
// ----------------

// Allowable level 16 images whose projection to level 8 give rise to QTT 8-congs
level16grps := [<grp[1], GL2Lift(grp[2], 16)> : grp in rzbdata | (GL2Level(grp[2]) le 16) and (grp[1] in all_QTT8_congs)]; 
level16grps := level16grps cat [<grp[1], GL2Project(grp[2], 16)> : grp in rzbdata | (GL2Level(grp[2]) gt 16) and (grp[1] in all_QTT8_congs)]; 

// Load level 16 data from Table 4.1
QTT163_congs := eval Read("table3-1data/QTT16-3congruencedata.txt");
all_QTT16_congs := {};

// These computations will take a while so we'll output as we go
printf "There are %o groups which work mod 8\n", #all_QTT8_congs;

// It takes a while to do G(ns 4+) so we'll do it seperately
assert level16grps[1][1] eq "X7";

printf "Working with group %o\n", level16grps[1][1];
// The quadratic twist must be by -Delta from the level 16 calculations. 
// Thus H is the full preimage of G(ns4) at level 16 so the following works

H := GL2Lift(GL2NonsplitCartan(4), 16);
H := GL2MinimizeGenerators(H);
assert #{g : g in Centralizer(GL(2, Integers(16)), H) | 2*Trace(g) eq 0} eq 0;
printf "-------------\n";

// For the rest of the groups
for grp in level16grps[2..#level16grps] do

    printf "Working with group %o\n", grp[1];

    flag, H_and_g := GL2GivesQTTCongruence(grp[2]);     
    for data in H_and_g do
        for power in {Integers()!Determinant(g) mod 8 : g in data[2]} do

            printf "Wow, it works mod 16 \n";

            assert power eq 3;
            assert exists(HH){HH : HH in QTT163_congs | IsInducedCongruence(HH[2], HH[3], grp[2], data[1])};

            Include(~all_QTT16_congs, grp[1]); 
        end for;
    end for;
    printf "-------------\n";
end for;


// ----------------
// LEVEL 32
// ----------------

assert all_QTT16_congs eq {"X55", "X441"}; // this is G(ns8+) and G(ns16+)

// ruling out the first case is easy, the quadratic twist must be by -Delta from
// the level 16 calculations. Thus H is the full preimage of G(ns8) at level 32

H := GL2Lift(GL2NonsplitCartan(8), 32);
H := GL2MinimizeGenerators(H);
assert #{g : g in Centralizer(GL(2, Integers(32)), H) | 2*Trace(g) eq 0} eq 0;
