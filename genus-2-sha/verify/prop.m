AttachSpec("../spec");

proofs := eval Read("../data/proof-data.m");
table := eval Read("../data/conductor-sorted-table.m");

//////////////////////////////////////////////////
// First prove the claim up of congruence quadratic twists
_<t> := PolynomialRing(Rationals());
j_X1 := (t^2 + t + 1)^(3)*(t^6 - 5*t^5 - 10*t^4 + 15*t^3 + 30*t^2 + 11*t + 1)^(3)/((t)^(7)*(t + 1)^(7)*(t^3 - 5*t^2 - 8*t - 1));

/* for proof in [p : p in proofs | p[3] ne []] do */
/*   // unpack genus 2 curve */
/*   APQ := proof[1]; */
/*   C := HyperellipticCurve(Genus2Curve(APQ)); */
/*   Kum := KummerSurface(Jacobian(C)); */
/*   f := DefiningPolynomial(Kum); */

/*   // unpack elliptic curve */
/*   E := EllipticCurve(proof[2]); */

/*   // unpack the fields K and L */
/*   OK := proof[3]; */
/*   K := NumberField(Polynomial(OK[2])); */
/*   OK := Order([Roots(Polynomial(x), K)[1][1] : x in OK]); */
  
/*   // unpack the torsion of Jac(C) into min polys of the torsion point coordinates */
/*   tors_J := proof[4]; */
/*   tors_J := [[OK!x : x in xi] : xi in tors_J]; */
/*   tors_J := [Polynomial(g) : g in tors_J]; */

/*   // the field Q(Kum[p]) */
/*   L := NumberField(tors_J[3]); */
/*   Kum_L := BaseChange(Kum, L); */

/*   // unpack the torsion of Jac(C) */
/*   xi := [[r[1] : r in Roots(t, L)] : t in tors_J]; */
/*   assert exists(P){Kum_L![1,xi1,xi2,xi3] : xi1 in xi[1], xi2 in xi[2], xi3 in xi[3] | Evaluate(f, [1,xi1,xi2,xi3]) eq 0}; */
/*   assert 7*P eq Kum_L![0,0,0,1]; */

/*   // unpack the torsion of E */
/*   tors_E := proof[5]; */
/*   tors_E := Polynomial([OK!x : x in tors_E]); */
/*   tors_E := Roots(tors_E, L)[1][1]; */
/*   assert Evaluate(j_X1, tors_E) eq jInvariant(E); */
/* end for; */


//////////////////////////////////////////////////
// Now prove the claim completely by checking that
// (1) Each pair has a proof of up to quad. twists
// (2) Traces of Frobenius pin down a unique quad. twist

for ex in table do
  APQ := ex[1];
  C := HyperellipticCurve(Polynomial(ex[2]));
  E := EllipticCurve(ex[3]);

  assert IsQuadraticTwist(C, HyperellipticCurve(Genus2Curve(APQ)));

  d := CorrectQuadraticTwist(C, E);
  assert d eq 1;
end for;


//////////////////////////////////////////////////
// Finally prove the claim proper about Sha[7]

fails := {@ @};

for ex in table do
  C := HyperellipticCurve(Polynomial(ex[2]));
  E := EllipticCurve(ex[3]);

  d := ex[4];
  C := QuadraticTwist(C, d);
  E := QuadraticTwist(E, d);
  C := ReducedMinimalWeierstrassModel(C);

  // check rank discrepency
  assert Rank(E) ge 2;
  _,ub := RankBounds(Jacobian(C));
  assert ub eq 0;

  // check bad primes
  assert GF(7)!(&*BadPrimes(C) * &*BadPrimes(E)) ne 0;

  // check Tamagawa numbers of E
  assert GF(7)!(&*TamagawaNumbers(E)) ne 0;

  // check Tamagawa numbers of C in Magma
  for p in BadPrimes(C) do
    try
      assert GF(7)!#ComponentGroup(RegularModel(C, p)) ne 0;
    catch e
      p, C;
      Include(~fails, <ex,p>);
      
      PrintFile("fails.m", fails : Overwrite:=true);

    end try;
  end for;  
end for;
