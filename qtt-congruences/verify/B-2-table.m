/*
    Proving that the genera listed in Table B.2 are correct.
*/

load "../functions/qttcong.m";

//------------
// SETUP
//------------
G3 := GL(2, Integers(3));
G4 := GL(2, Integers(4));
G5 := GL(2, Integers(5));
G7 := GL(2, Integers(7));
G9 := GL(2, Integers(9));
G11 := GL(2, Integers(11));

Gns3 := GL2NonsplitCartan(3);
Gs3 := GL2SplitCartan(3);
Gns5 := GL2NonsplitCartan(5);
Gs5 := GL2SplitCartan(5);
Gns7 := GL2NonsplitCartan(7);
Gs7 := GL2SplitCartan(7);
Gns9 := GL2NonsplitCartan(9);
Gns11 := GL2NonsplitCartan(11);


columns :=[
            <Gns3, Normaliser(G3, Gns3)>,
            <Gs3, Normaliser(G3, Gs3)>,
            <Gns5, Normaliser(G5, Gns5)>,
            <Gs5, Normaliser(G5, Gs5)>,
            <Gns7, Normaliser(G7, Gns7)>,
            <Gs7, Normaliser(G7, Gs7)>,
            <Gns9, Normaliser(G9, Gns9)>,
            <Gns11, Normaliser(G11, Gns11)>
];

QTT41_congs := eval Read("table3-1data/QTT4-1congruencedata.txt"); 
QTT43_congs := eval Read("table3-1data/QTT4-3congruencedata.txt");
QTT81_congs := eval Read("table3-1data/QTT8-1congruencedata.txt"); 
QTT83_congs := eval Read("table3-1data/QTT8-3congruencedata.txt");
QTT85_congs := eval Read("table3-1data/QTT8-5congruencedata.txt"); 
QTT87_congs := eval Read("table3-1data/QTT8-7congruencedata.txt");
QTT163_congs := eval Read("table3-1data/QTT16-3congruencedata.txt"); 

rows := [
        <QTT41_congs[1][3], QTT41_congs[1][2]>,
        <QTT41_congs[2][3], QTT41_congs[2][2]>,
        <QTT43_congs[1][3], QTT43_congs[1][2]>,
        <QTT43_congs[2][3], QTT43_congs[2][2]>,  
        <QTT43_congs[3][3], QTT43_congs[3][2]>,
        <QTT81_congs[1][3], QTT81_congs[1][2]>,
        <QTT83_congs[1][3], QTT83_congs[1][2]>,
        <QTT83_congs[2][3], QTT83_congs[2][2]>,  
        <QTT85_congs[1][3], QTT85_congs[1][2]>,
        <QTT87_congs[1][3], QTT87_congs[1][2]>,
        <QTT87_congs[2][3], QTT87_congs[2][2]>,
        <QTT163_congs[1][3], QTT163_congs[1][2]>
];



tableA_2 := [[1,1,5,7,13,17,19,41],
          [1,3,7,9,15,21,21,45],
          [0,0,0,1,1,1,1,4],
          [2,3,7,10,16,21,22,46],
          [0,1,5,8,12,17,16,40],
          [7,-1,-1,-1,-1,-1,-1,-1],
          [0,1,3,6,8,11,10,-1],
          [5,13,-1,-1,-1,-1,-1,-1],
          [7,-1,-1,-1,-1,-1,-1,-1],
          [2,5,13,20,30,41,40,-1],
          [5,13,-1,-1,-1,-1,-1,-1],
          [3,9,-1,-1,-1,-1,-1,-1]];


//------------
// CREATE THE TABLE
//------------

table := [];
for HHp_r in rows[1..5] do 
    curr_row := [];
    for HHp_c in columns do 
        Hdelta := GL2ModDeltaGroup(HHp_r[1], HHp_r[2], HHp_c[1], HHp_c[2]);
        Append(~curr_row, GL2Genus(Hdelta));
    end for;
    Append(~table, curr_row);
end for;

// We need to know what to ignore at level 8. Note that the image of
// G(ns{2^{k+1}}^+) at level 2^k is contained in G(ns{2^k}^+) and the image of 
// G(s4^+) at level 2 is contained in G(ns2^+).

assert IsConjugate(G4, GL2Project(QTT81_congs[1][2], 4), QTT41_congs[2][2]); // Row 6: pi_4(G_97) = G_11 <- row 2
assert IsConjugate(G4, GL2Project(QTT83_congs[2][2], 4), QTT43_congs[3][2]); // Row 8: pi_4(G_90) = G_11 <- row 5
assert IsConjugate(G4, GL2Project(QTT85_congs[1][2], 4), QTT41_congs[2][2]); // Row 9: pi_4(G_97) = G_11 <- row 2
assert IsConjugate(G4, GL2Project(QTT87_congs[2][2], 4), QTT43_congs[3][2]); // Row 11: pi_4(G_90) = G_11 <- row 5

assert [i : i in [1..8] | table[2][i] le 1] eq [1];
assert [i : i in [1..8] | table[5][i] le 1] eq [1,2];
assert [i : i in [1..8] | table[3][i] le 1] eq [1..7];

// So now we can do the rest of the table ignoring the appropriate columns.

// G_97
HHp_r := rows[6];
curr_row := [];
for HHp_c in columns[1..1] do 
    Hdelta := GL2ModDeltaGroup(HHp_r[1], HHp_r[2], HHp_c[1], HHp_c[2]);
    Append(~curr_row, GL2Genus(Hdelta));
end for;
curr_row := curr_row cat [-1 : i in [1..7]];
Append(~table, curr_row);

// G(ns4+)
HHp_r := rows[7];
curr_row := [];
for HHp_c in columns[1..7] do 
    Hdelta := GL2ModDeltaGroup(HHp_r[1], HHp_r[2], HHp_c[1], HHp_c[2]);
    Append(~curr_row, GL2Genus(Hdelta));
end for;
curr_row := curr_row cat [-1];
Append(~table, curr_row);

// G_90
HHp_r := rows[8];
curr_row := [];
for HHp_c in columns[1..2] do 
    Hdelta := GL2ModDeltaGroup(HHp_r[1], HHp_r[2], HHp_c[1], HHp_c[2]);
    Append(~curr_row, GL2Genus(Hdelta));
end for;
curr_row := curr_row cat [-1 : i in [1..6]];
Append(~table, curr_row);

// G_63
HHp_r := rows[9];
curr_row := [];
for HHp_c in columns[1..1] do 
    Hdelta := GL2ModDeltaGroup(HHp_r[1], HHp_r[2], HHp_c[1], HHp_c[2]);
    Append(~curr_row, GL2Genus(Hdelta));
end for;
curr_row := curr_row cat [-1 : i in [1..7]];
Append(~table, curr_row);

// G(s4+)
HHp_r := rows[10];
curr_row := [];
for HHp_c in columns[1..7] do 
    Hdelta := GL2ModDeltaGroup(HHp_r[1], HHp_r[2], HHp_c[1], HHp_c[2]);
    Append(~curr_row, GL2Genus(Hdelta));
end for;
curr_row := curr_row cat [-1];
Append(~table, curr_row);

// G_90
HHp_r := rows[11];
curr_row := [];
for HHp_c in columns[1..2] do 
    Hdelta := GL2ModDeltaGroup(HHp_r[1], HHp_r[2], HHp_c[1], HHp_c[2]);
    Append(~curr_row, GL2Genus(Hdelta));
end for;
curr_row := curr_row cat [-1 : i in [1..6]];
Append(~table, curr_row);

// finally at level 16 we want to ignore the last 6 columns
assert [i : i in [1..8] | table[7][i] in {0,1}] eq [1,2]; // the G(ns4+) row

HHp_r := rows[12];
curr_row := [];
for HHp_c in columns[1..2] do 
    Hdelta := GL2ModDeltaGroup(HHp_r[1], HHp_r[2], HHp_c[1], HHp_c[2]);
    Append(~curr_row, GL2Genus(Hdelta));
end for;
curr_row := curr_row cat [-1 : i in [1..6]];
Append(~table, curr_row);

//------------
// CHECK TABLES ARE EQUAL
//------------

assert table eq tableA_2;
