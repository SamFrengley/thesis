/*
    Checking that the groups in Table 4.1 all have H (and H+) conjugate to its
    transpose.
*/

QQ := Rationals();
Attach("../RSZ21/gl2.m");

// ----------------
// LEVEL 4
// ----------------

// Load level 4 data from Table 4.1
G := GL(2, Integers(4));
QTT41_congs := eval Read("table3-1data/QTT4-1congruencedata.txt"); 
QTT43_congs := eval Read("table3-1data/QTT4-3congruencedata.txt");

for data in QTT41_congs do
    assert IsConjugate(G, data[2], GL2Transpose(data[2]));
    assert IsConjugate(G, data[3], GL2Transpose(data[3]));
end for;

for data in QTT43_congs do
    assert IsConjugate(G, data[2], GL2Transpose(data[2]));
    assert IsConjugate(G, data[3], GL2Transpose(data[3]));
end for;

// ----------------
// LEVEL 8
// ----------------

// Load level 8 data from Table 4.1
G := GL(2, Integers(8));
QTT81_congs := eval Read("table3-1data/QTT8-1congruencedata.txt");
QTT83_congs := eval Read("table3-1data/QTT8-3congruencedata.txt");
QTT85_congs := eval Read("table3-1data/QTT8-5congruencedata.txt"); 
QTT87_congs := eval Read("table3-1data/QTT8-7congruencedata.txt");

for data in QTT81_congs do
    assert IsConjugate(G, data[2], GL2Transpose(data[2]));
    assert IsConjugate(G, data[3], GL2Transpose(data[3]));
end for;

for data in QTT83_congs do
    assert IsConjugate(G, data[2], GL2Transpose(data[2]));
    assert IsConjugate(G, data[3], GL2Transpose(data[3]));
end for;

for data in QTT85_congs do
    assert IsConjugate(G, data[2], GL2Transpose(data[2]));
    assert IsConjugate(G, data[3], GL2Transpose(data[3]));
end for;

for data in QTT87_congs do
    assert IsConjugate(G, data[2], GL2Transpose(data[2]));
    assert IsConjugate(G, data[3], GL2Transpose(data[3]));
end for;

// ----------------
// LEVEL 16
// ----------------

// Load level 16 data from Table 4.1
G := GL(2, Integers(16));
QTT163_congs := eval Read("table3-1data/QTT16-3congruencedata.txt");

for data in QTT163_congs do
    assert IsConjugate(G, data[2], GL2Transpose(data[2]));
    assert IsConjugate(G, data[3], GL2Transpose(data[3]));
end for;

// ----------------
// LEVEL 32
// ----------------
G := GL(2, Integers(32));
QTT323_congs := eval Read("table3-1data/QTT32-3congruencedata.txt");

for data in QTT323_congs do
    assert IsConjugate(G, data[2], GL2Transpose(data[2]));
    assert IsConjugate(G, data[3], GL2Transpose(data[3]));
end for;
