/*
    Proving that the genera listed in Table B.1 are correct.
*/

load "../functions/qttcong.m";

//------------
// SETUP
//------------
G3 := GL(2, Integers(3));
G5 := GL(2, Integers(5));
G7 := GL(2, Integers(7));
G11 := GL(2, Integers(11));

Gns3 := GL2NonsplitCartan(3);
Gs3 := GL2SplitCartan(3);
Gns5 := GL2NonsplitCartan(5);
Gs5 := GL2SplitCartan(5);
Gns7 := GL2NonsplitCartan(7);
Gs7 := GL2SplitCartan(7);
Gns11 := GL2NonsplitCartan(11);

columns := [
            <Gns5, Normaliser(G5, Gns5)>,
            <Gs5, Normaliser(G5, Gs5)>,
            <Gns7, Normaliser(G7, Gns7)>,
            <Gs7, Normaliser(G7, Gs7)>,
            <Gns11, Normaliser(G11, Gns11)>
];

rows :=[
            <Gns3, Normaliser(G3, Gns3)>,
            <Gs3, Normaliser(G3, Gs3)>,
            <Gns5, Normaliser(G5, Gns5)>,
            <Gs5, Normaliser(G5, Gs5)>,
            <Gns7, Normaliser(G7, Gns7)>,
            <Gs7, Normaliser(G7, Gs7)>
];

tableA_1 := [[2,3,4,7,17],
          [5,8,12,17,40],
          [-1,-1,26,35,76],
          [-1,-1,40,55,117],
          [-1,-1,-1,-1,166],
          [-1,-1,-1,-1,225]];


//------------
// CREATE THE TABLE
//------------

table := [];
for HHp_r in rows do 
    curr_row := [];
    for HHp_c in columns do 
        p1 := GL2Level(HHp_r[1]); p2 := GL2Level(HHp_c[2]);
        if p1 ge p2 then
            Append(~curr_row, -1);
        else
            Hdelta := GL2ModDeltaGroup(HHp_r[1], HHp_r[2], HHp_c[1], HHp_c[2]);
            Append(~curr_row, GL2Genus(Hdelta));
        end if;
    end for;
    Append(~table, curr_row);
end for;


assert table eq tableA_1;
