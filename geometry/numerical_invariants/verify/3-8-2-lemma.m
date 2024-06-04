AttachSpec("../spec");

function Possible_r(N)
  all := [r : r in [1..N] | GCD(N, r) eq 1];
  mod_sq := [];
  for r in all do
    if not exists{s : s in mod_sq | IsSquare(Integers(N)!r*s)} then
      Append(~mod_sq, r);
    end if;
  end for;
  return mod_sq;
end function;

//////////////////////////////////////////////////
// Proof of Lemma 3.8.2
//////////////////////////////////////////////////
for N in [23] cat [25..105] do
  for r in Possible_r(N) do
    assert KWoSquared(N, r) gt 0;
  end for;
end for;

assert KWoSquared(21,1) gt 0;
assert KWoSquared(21,10) gt 0;
assert KWoSquared(24,13) gt 0;
assert KWoSquared(24,19) gt 0;
