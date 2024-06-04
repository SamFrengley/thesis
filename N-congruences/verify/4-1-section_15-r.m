Attach("../ZNr-equations.m");
QQ := Rationals();


function CorrectMap(m, ff)
// This returns the map from X_0(m) (in the small modular curve database) to the copy of 
// X_0^+(m) on our model for Z(15, m). This will be used to check the curves we have 
// do map to Hecke correspondences
  if m eq 8 then //blown down
    x0 := Explode(ff);
    pp := [(x0^2 + 11*x0 + 32)/x0];

  elif m eq 11 then //blown down
    x0,x1 := Explode(ff);
    pp := [1/(x0 - 16)*x1 + (3*x0 + 13)/(x0 - 16)];
    
  elif m eq 17 then
    x0,x1 := Explode(ff);
    pp := [
      (-1/2*x0 + 7/2)/(x0 + 1/2*x1 + 7/2),
      1
    ];

  elif m eq 23 then
    x0,x1 := Explode(ff);
    pp := [
      1/x0,
      0
    ];

  elif m eq 26 then
    x0,x1 := Explode(ff);
    pp := [
      0,
      -1/x0
    ];

  elif m eq 29 then
    x0,x1 := Explode(ff);
    pp := [
      -1/x0,
      x0
    ];

  elif m eq 32 then
    x0,x1 := Explode(ff);
    pp := [
      (-x0^2 + 2*x0)/(x0*x1 + 2*x0 + x1),
      (2*x0 + x1)/(x0 + x1 + 2)
    ];

  elif m eq 38 then
    x0,x1,x2 := Explode(ff);
    pp := [
      (x0 - 1)/(x0 - x1),
      (x2 + 1)/(x0 - x1 + x2 + 1)
    ];

  elif m eq 41 then
    x0,x1 := Explode(ff);
    pp := [
      1/x0,
      -1/x0
    ];

  elif m eq 44 then
    x0,x1,x2 := Explode(ff);
    pp := [
      (x0 + x2)/(x1 - x2 -1),
      1/(x0 + x2)
    ];      
    
  elif m eq 47 then
    x0,x1 := Explode(ff);
    pp := [
      -1/(x0 - 1),
      x0/(x0 - 1)
    ];

  elif m eq 53 then
    x0,x1,x2 := Explode(ff);
    pp := [
      -x2/(x0 - x2 + 1),
      (x2 - 1)/x0
    ];

  elif m eq 56 then
    x0,x1 := Explode(ff);
    pp := [
      (-x0^2 + 2*x0*x1 + x1^2 - 2)/(x0^2*x1 - x0*x1^2 - x0 + x1),
      (x0^2 - 2*x0*x1 - x1^2 + 2)/(x0^3 - x0^2 - x0*x1^2 + x1^2)
    ];

  elif m eq 59 then
    x0,x1 := Explode(ff);
    pp := [
      -1/x0,
      (x0 - 1)
    ];

  elif m eq 71 then
    x0,x1 := Explode(ff);
    pp := [
      1/x0,
      -1/(x0 + 1)
    ];
    
  end if;
  return pp;
end function;


//////////////////////////////////////////////////////////////////////
A3<u,v,z> := AffineSpace(QQ, 3);
K<x0,x1,x2> := FunctionField(QQ, 3);
A1 := AffineSpace(QQ, 1);
A2<a,b> := AffineSpace(QQ, 2);
Kt<t> := FunctionField(QQ);
P<e> := PowerSeriesRing(Kt);

/**************************************************
  Z(15,2)
*/
ZZ := Surface(A3, ZNrEquations(15,2 : vars:=[u,v,z]));
print "--------------------\nThe case Z(15,2)\n";

////////////////////
// FIRST DEAL WITH THE BLOWN DOWN CASES
print "The blown down cases:\n";
f := -ZNrEquations(15,2 : vars:=[x0/x2,x1/x2,0])[1];
f := x2^14 * f;

// X_0(8)
printf "The case when m=8\n";
bl_up_1 := Evaluate(f, [e,-e+e^2*t,1]);
assert bl_up_1 + BigO(e^5) eq (t^2 - 22*t - 7)*e^4 + BigO(e^5);
bl_up_4 := Evaluate(f, [e^2*t,1,e]);
assert bl_up_4 + BigO(e^5) eq (t^2 - 24*t + 16)*e^4 + BigO(e^5);

X0m := SmallModularCurve(8);
X0m_KK := BaseChange(X0m, Kt);
j1 := jNInvariant(X0m_KK![t,1], 8);
j2 := jInvariant(X0m_KK![t,1], 8);

//case 1
pp := CorrectMap(8, [t]);
jj,jj_1728 := WNrModuli(15,2 : vars:=[e,-e+e^2*t]);
jj := Evaluate(jj,0); jj_1728 := Evaluate(jj_1728,0); // On the exceptional divisor
jj := Evaluate(jj, pp[1]); jj_1728 := Evaluate(jj_1728, pp[1]);
assert jj eq j1*j2;
assert jj_1728 eq (j1 - 1728)*(j2 - 1728);

//case 4
pp := [(t^2 + 12*t + 32)/t];
jj,jj_1728 := WNrModuli(15,2 : vars:=[e*t,1/e]);
jj := Evaluate(jj,0); jj_1728 := Evaluate(jj_1728,0); // On the exceptional divisor
jj := Evaluate(jj, pp[1]); jj_1728 := Evaluate(jj_1728, pp[1]);
assert jj eq j1*j2;
assert jj_1728 eq (j1 - 1728)*(j2 - 1728);

// X_0(17)
printf "The case when m=17\n";
bl_up_1 := Evaluate(f, [1,e + e^2 + e^3*t,e]);
assert bl_up_1 + BigO(e^15) eq -(16*t^4 - 92*t^3 + 207*t^2 - 208*t + 76)*e^14 + BigO(e^15);

X0m := SmallModularCurve(17);
KK<x0,x1> := FunctionField(X0m);

pp := [-1/(x0^2 + 4*x0 + 8)*x1 + (x0^2 + 5*x0 + 14)/(x0^2 + 4*x0 + 8)];

jj,jj_1728 := WNrModuli(15,2 : vars:=[1/e,1 + e + e^2*t]);
jj := Evaluate(jj,0); jj_1728 := Evaluate(jj_1728,0); // On the exceptional divisor
jj := Evaluate(jj, pp[1]); jj_1728 := Evaluate(jj_1728, pp[1]);

X0m_KK := BaseChange(X0m, KK);
j1 := jNInvariant(X0m_KK![x0,x1,1], 17);
j2 := jInvariant(X0m_KK![x0,x1,1], 17);
assert jj eq j1*j2;
assert jj_1728 eq (j1 - 1728)*(j2 - 1728);

////////////////////
// NOW THE NON-BLOWN DOWN
print "\nAnd the remaining cases:\n";
table := eval Read("../Z15-r/curvesonZ15-r/modularcurves/15-2.m");

for ex in table do
  m, X, _ := Explode(ex);
  printf "The case when m=%o\n", m;
  X := Curve(ZZ, [X]);
  X0m := SmallModularCurve(m);
  g := Genus(X0m);
  assert Genus(X) eq g;

  if g eq 0 then
    KK := FunctionField(QQ);
    coords := [KK.1];
    gen_pt := [KK.1, 1];
  else
    KK := FunctionField(X0m);
    PP := AmbientSpace(X0m);
    gr := Gradings(PP)[1];
    N := #gr;
    assert gr[N] eq 1;
    gr := [Integers()!(g/gr[N]) : g in gr];
    gen_pt := [KK!(X0m.i/(X0m.N)^gr[i]) : i in [1..N]];
    coords := gen_pt[1..N-1];
    X0m := BaseChange(X0m, KK);    
  end if;

  pp := CorrectMap(m, coords);

  jj, jj_1728 := WNrModuli(15,2 : vars:=pp);
  j1 := jNInvariant(X0m!gen_pt, m);
  j2 := jInvariant(X0m!gen_pt, m);

  assert jj eq j1*j2;
  assert jj_1728 eq (j1 - 1728)*(j2 - 1728);
end for;


/**************************************************
  Z(15,11)
*/
A3<u,v,z> := AffineSpace(QQ, 3);
K<x0,x1,x2> := FunctionField(QQ, 3);
A1 := AffineSpace(QQ, 1);
A2<a,b> := AffineSpace(QQ, 2);
Kt<t> := FunctionField(QQ);
P<e> := PowerSeriesRing(Kt);

print "--------------------\nThe case Z(15,11)\n";
ZZ := Surface(A3, ZNrEquations(15,11 : vars:=[u,v,z]));

////////////////////
// FIRST DEAL WITH THE BLOWN DOWN CASES
print "The blown down cases:\n";
f := -ZNrEquations(15,11 : vars:=[x0/x2,x1/x2,0])[1];
f := x2^12 * f;

// X_0(11)
printf "The case when m=11\n";
bl_up := Evaluate(f, [-1 + e^2*t, 1, e]);
assert bl_up + BigO(e^3) eq (t^4 - 12*t^3 - 40*t^2 - 28*t - 8)*e^2 + BigO(e^3);

X0m := SmallModularCurve(11);
KK<x0,x1> := FunctionField(X0m);

// Check that this curve is isomorphic to the exceptional divisor
pp := [
  1/(x0 - 16)*x1 + (3*x0 + 13)/(x0 - 16),
  121/(x0^2 - 32*x0 + 256)*x1 + (-x0^3 + 48*x0^2 - 42*x0 - 139)/(x0^2 - 32*x0 + 256)
];
F_11 := Curve(A2, b^2 - (a^4 - 12*a^3 - 40*a^2 - 28*a - 8));
psi := map<X0m -> F_11 | pp>;
assert IsInvertible(psi);

// Check the moduli interpretation
jj,jj_1728 := WNrModuli(15,11 : vars:=[(-1 + e^2*t)/e, 1/e]);
jj := Evaluate(jj,0); jj_1728 := Evaluate(jj_1728,0); // On the exceptional divisor
jj := Evaluate(jj, pp[1]); jj_1728 := Evaluate(jj_1728, pp[1]);

X0m_KK := BaseChange(X0m, KK);
j1 := jNInvariant(X0m_KK![x0,x1,1], 11);
j2 := jInvariant(X0m_KK![x0,x1,1], 11);

assert j1*j2 eq jj;
assert (j1-1728)*(j2-1728) eq jj_1728;

////////////////////
// NOW THE NON-BLOWN DOWN
print "\nAnd the remaining cases:\n";
table := eval Read("../Z15-r/curvesonZ15-r/modularcurves/15-11.m");

for ex in table do
  m, X, _ := Explode(ex);
  printf "The case when m=%o\n", m;
  X := Curve(ZZ, [X]);
  X0m := SmallModularCurve(m);
  g := Genus(X0m);
  assert Genus(X) eq g;

  if g eq 0 then
    KK := FunctionField(QQ);
    coords := [KK.1];
    gen_pt := [KK.1, 1];
  else
    KK := FunctionField(X0m);
    PP := AmbientSpace(X0m);
    gr := Gradings(PP)[1];
    N := #gr;
    assert gr[N] eq 1;
    gr := [Integers()!(g/gr[N]) : g in gr];
    gen_pt := [KK!(X0m.i/(X0m.N)^gr[i]) : i in [1..N]];
    coords := gen_pt[1..N-1];
    X0m := BaseChange(X0m, KK);    
  end if;

  pp := CorrectMap(m, coords);

  jj, jj_1728 := WNrModuli(15,11 : vars:=pp);
  j1 := jNInvariant(X0m!gen_pt, m);
  j2 := jInvariant(X0m!gen_pt, m);

  assert jj eq j1*j2;
  assert jj_1728 eq (j1 - 1728)*(j2 - 1728);
end for;




