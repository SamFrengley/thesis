Attach("../../../ZNr-equations.m");
prime := NextPrime(10^10);
F := GF(prime);
A2<u,v> := AffineSpace(F, 2);
P<[x]> := PolynomialRing(Rationals(), 8);
jj, jj_1728 := WNrModuli(14,3 : vars:=[u,v]);

crvs := eval Read("14-3.m");

k := 1;

for i in [6..#crvs] do
  crv := crvs[i];
  m := crv[1]; m;
  C := Curve(A2, crv[k+1]);
  j_map := map<C -> A2 | [jj,jj_1728]>;
  
  X := SmallModularCurve(m);
  KX<[xx]> := FunctionField(X);
  //Image(map<X -> AffineSpace(Rationals(), 2) | [-xx[3]/(xx[1] - xx[3] + 1), (xx[3] - 1)/(xx[1])]>);
  X_K := BaseChange(X, KX);
  
  PP := AmbientSpace(X);
  gr := Gradings(PP)[1];
  N := #gr;
  assert gr[N] eq 1;
  gr := [Integers()!(g/gr[N]) : g in gr];
  gen_pt := [KX!(X.i/(X.N)^gr[i]) : i in [1..N]];
  
  j1 := jInvariant(X_K!gen_pt, m);
  j2 := jNInvariant(X_K!gen_pt, m);

  /////////
  X := ChangeRing(X, F);
  _<[xx]> := FunctionField(X);
  j1 := eval Sprint(j1);
  j2 := eval Sprint(j2);

  pts := GeneratePoints(X, 100);
  pts := [X!pt : pt in pts];

  inits := [];
  images := [];

  bad := BasePoints(j_map);
  for pt in pts do
    W1_pt := A2![j1(pt)*j2(pt), (j1(pt)-1728)*(j2(pt) - 1728)];
    ims := [p : p in Points((W1_pt @@ j_map)) | not p in bad];
    if #ims eq 1 then
      Append(~images, Eltseq(ims[1]));
      Append(~inits, Eltseq(pt));
    end if;
  end for;

  found := false;
  d := 0;
  while not found do
    d := d + 1;
    m1 := InterpolateMapToA1(inits, [im[1] : im in images], d);
    if #m1 ge 1 then
      assert exists(m1){x : x in m1 | x[2] ne 0};
      found := true;
    end if;
  end while;

  found := false;
  d := 0;
  while not found do
    d := d + 1;
    m2 := InterpolateMapToA1(inits, [im[2] : im in images], d);
    if #m2 ge 1 then
      assert exists(m2){x : x in m2 | x[2] ne 0};
      found := true;
    end if;
  end while;

  m := [m1[1]/m1[2], m2[1]/m2[2]];
  m;
  [x : x in m1]; [x : x in m2];
  Image(map<X -> A2 | m>);
  
end for;
