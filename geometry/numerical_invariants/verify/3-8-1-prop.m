AttachSpec("../spec");

//////////////////////////////////////////////////
// Proof of Proposition 3.8.1
//////////////////////////////////////////////////

assert barCinf_barKW(6,5) le -2;
assert barCinf_barKW(7,6) le -2;
assert barCinf_barKW(8,3) le -2;
assert barCinf_barKW(8,5) le -2;
assert barCinf_barKW(8,7) le -2;

for N in [9,10] do
  for r in [r : r in [1..N] | GCD(N,r) eq 1] do
    assert barCinf_barKW(N,r) le -2;
  end for;
end for;

assert barCinf_barKW(12,5) le -1;
assert barCinf_barKW(12,7) le -1;
assert barCinf_barKW(12,11) le -1;
assert barCinf_barKW(13,1) le -1;
assert barCinf_barKW(13,2) le -1;
assert barCinf_barKW(14,13) le -1;
assert barCinf_barKW(15,2) le -1;
assert barCinf_barKW(15,11) le -1;
assert barCinf_barKW(16,3) le -1;
assert barCinf_barKW(16,7) le -1;
assert barCinf_barKW(18,5) le -1;
assert barCinf_barKW(20,11) le -1;
