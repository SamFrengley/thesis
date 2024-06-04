Attach("../ZNr-equations.m");
QQ := Rationals();

function CorrectMap(m, ff)
// This returns the map from X_0(m) (in the small modular curve database) to the copy of 
// X_0^+(m) on our model for Z(14, m). This will be used to check the curves we have 
// do map to Hecke correspondences
  if m eq 25 then
    x0 := Explode(ff);
    pp := [
      (x0^2 + 3*x0 + 5)/x0,
      (x0^2 + 3*x0 + 5)/x0 + 1
    ];

  elif m eq 27 then
    x0,x1 := Explode(ff);
    pp := [
      2*(-3*x0^2 - 7*x0*x1 + 4*x0 - 2*x1^2 - 6*x1 - 7)/(5*x0^2 + 5*x0*x1 - 8*x0 + x1^2 - 4*x1),
      (-4*x1 - 20)/(8*x0 + 4*x1 - 4)
    ];
    
  elif m eq 29 then
    x0,x1 := Explode(ff);
    pp := [
      -x0,
      (x0 - 1)/x0
    ];

  elif m eq 31 then
    x0,x1 := Explode(ff);
    pp := [
      (2*x0 + 2)/(-x0^2 + 2*x0 - 1),
      (-x0 - 1)/(x0 - 1)
    ];

  elif m eq 33 then
    x0,x1 := Explode(ff);
    pp := [
      (24*x0^4 -20*x0^3 - 2*x0^2*x1 + 72*x0^2 - 2*x0*x1^2 + 4*x0*x1 - 36*x0 -
        4*x1^2 + 4*x1 + 32)/(2*x0^4 + x0^3*x1 + 3*x0^2*x1 + x0*x1^2 - x0*x1 - 2*x0),
      (-3*x0^2 + 3*x0 - x1 - 3)/(-x0^2 - x0 - x1 - 1)
    ];
    
  elif m eq 37 then 
    x0,x1 := Explode(ff);
    pp := [
      (2*x0^3 - x0*x1 - 1)/(6*x0^3 - 2*x0^2*x1 - 3*x0^2 + 3*x0 - x1^2),
      (x0^2 - x0)/(x1)
    ];

  elif m eq 39 then
    x0,x1 := Explode(ff);
    pp := [
      x0,
      x0^2/(x0 - 1)
    ];

  elif m eq 41 then
    x0,x1 := Explode(ff);
    pp := [
      (4*x0^2 - 2*x0 - 2)/(-x0^2 - x0),
      (x0 - 2)/(-x0)
    ];
    
  elif m eq 43 then
    x0,x1 := Explode(ff);
    pp := [
      (x0^2 + x0*x1 + x1^2 + 2)/(x1),
      (x0^2 + x0*x1 + 2)/(x1)
    ];

  elif m eq 45 then
    x0,x1 := Explode(ff);
    pp := [
      (4*x0^2 + 2*x0*x1 - 6*x0 - 2*x1 + 2)/(-x0^2 + x0*x1 + 2*x0 - x1^2 - 2*x1 - 4),
      (-x0^2 - 5*x0*x1 + 4*x0 - x1^2 - 4*x1 + 10)/(x0^2 + 3*x0*x1 + 2*x0 + x1^2 - 2*x1)
    ];

  elif m eq 47 then
    x0,x1 := Explode(ff);
    pp := [
      (-2*x0^3 - 2)/(x0^3 - 2*x0^2 + x0),
      (x0^2 - x0 + 2)/(x0^2 - x0)
    ];

  elif m eq 51 then
    x0,x1 := Explode(ff);
    pp := [
      (x0^2 - x0*x1 - x0 + 1)/(-2*x0^2 - x1),
      (x0^3 + x0^2 - x0*x1^2 - x0 - x1^2 + x1)/(2*x0^3 +
        x0^2*x1 - 3*x0^2 - x0*x1^2 + x0*x1 + x0 - 1)
    ];
    
  elif m eq 53 then
    x0,x1,x2 := Explode(ff);
    pp := [
      -1/(x2 - 1),
      (-x0 + x2)/x2
    ];
    
  elif m eq 55 then
    x0,x1 := Explode(ff);
    pp := [
      (2*x0^2*x1 + 10*x0*x1^2 - 4*x0 + 10*x1^2 + 8*x1 - 4)/(x0^2*x1 - 4*x0*x1^2 + 
        3*x0*x1 - 4*x1^2 + 7*x1),
      (x0^2 + 2*x0*x1 + 5*x0 + 2*x1 + 9)/(x0^2 - 4*x0*x1 + 3*x0 - 4*x1 + 7)
    ];
    
  elif m eq 57 then
    x0,x1 := Explode(ff);
    pp := [ 
      (-x0^2 - x0*x1^2 - x0*x1 + 2*x0 - x1^3 - x1^2 - 5*x1 - 1)/(x0^2 -
        x0*x1^2 - 2*x0*x1 + x0 - x1^2 - x1 - 2),
      (x0^2 + x0*x1^2 + 2*x0*x1 + 4*x0 + x1^2 + x1 + 4)/(x0^2 + x0*x1^2 +
        4*x0 + x1^2 + 4)
    ];

  elif m eq 59 then
    x0,x1 := Explode(ff);
    pp := [
      (-2*x0^2 - 2*x0 + 2)/(x0^3 - 2*x0^2 + 2*x0 - 1),
      (-x0 - 1)/(x0 - 1)
    ];

  elif m eq 61 then
    x0, x1,x2 := Explode(ff);
    pp := [
      (-4*x0*x1 + 4*x1^2 - 4*x1*x2 - 2)/(x0*x1 - x1^2 + x1*x2 - 2*x1),
      (x0 - x1 + x2 - 1)/(x0 + x1 + x2 - 1)
    ];
    
  elif m eq 71 then
    x0,x1 := Explode(ff);
    pp := [
      -x0,
      (x0^3 - x0)/(x0^3 + x0^2 - 1)
    ];
    
  elif m eq 81 then
    x0,x1,x2 := Explode(ff);
    pp := [
      1/(x0 - x1),
      (3*x0*x1 - x0*x2 - x0 + x1)/(x0^2 + x0*x1 + x1^2 - x1*x2)
    ];
    
  end if;

  return pp;
end function;

//////////////////////////////////////////////////////////////////////
A3<u,v,z> := AffineSpace(QQ, 3);
K<[x]> := FunctionField(QQ, 3);
A1 := AffineSpace(QQ, 1);
A2<a,b> := AffineSpace(QQ, 2);
Kt<t> := FunctionField(QQ);
P<e> := PowerSeriesRing(Kt);

/**************************************************
  Z(14,1)
*/
ZZ := Surface(A3, ZNrEquations(14,1 : vars:=[u,v,z]));
print "--------------------\nThe case Z(14,1)\n";

////////////////////
// FIRST DEAL WITH THE BLOWN DOWN CASES
print "The blown down cases:\n";
f := -ZNrEquations(14,1 : vars:=[x[1]/x[3],x[2]/x[3],0])[1];
f := x[3]^12 * f;

// X_0(9)
printf "The case when m=9\n";
bl_up := Evaluate(f, [e,-e + e^2*t,1]);
assert bl_up + BigO(e^5) eq (t^2 - 16*t - 44)*e^4 + BigO(e^5);

// Parametrise the exceptional divisor (where e=0)
pp := [
  (t^2 + 8*t + 27)/t,
  (t^2 - 27)/t
];
F_9 := Curve(A2, b^2 - (a^2 - 16*a - 44));
psi := map<A1 -> F_9 | [Evaluate(p, A1.1) : p in pp]>;
assert IsInvertible(psi);

// Check the moduli interpretation
jj,jj_1728 := WNrModuli(14,1 : vars:=[e,-e + e^2*t]);
jj := Evaluate(jj,0); jj_1728 := Evaluate(jj_1728,0); // On the exceptional divisor
jj := Evaluate(jj, pp[1]); jj_1728 := Evaluate(jj_1728, pp[1]);

X0m := SmallModularCurve(9);
X0m := BaseChange(X0m, Kt);
j1 := jNInvariant(X0m![t,1], 9);
j2 := jInvariant(X0m![t,1], 9);

assert j1*j2 eq jj;
assert (j1-1728)*(j2-1728) eq jj_1728;

// X_0(11)
printf "The case when m=11\n";
bl_up := Evaluate(f, [1 + e^2*t, 1, e]);
assert bl_up + BigO(e^9) eq (t^4 + 12*t^3 - 40*t^2 + 28*t - 8)*e^8 + BigO(e^9);

X0m := SmallModularCurve(11);
KK<x0,x1> := FunctionField(X0m);

// Check that this curve is isomorphic to the exceptional divisor
pp := [
  -1/(x0 - 16)*x1 + (-3*x0 - 13)/(x0 - 16),
  121/(x0^2 - 32*x0 + 256)*x1 + (-x0^3 + 48*x0^2 - 42*x0 - 139)/(x0^2 - 32*x0 + 256)
];
F_11 := Curve(A2, b^2 - (a^4 + 12*a^3 - 40*a^2 + 28*a - 8));
psi := map<X0m -> F_11 | pp>;
assert IsInvertible(psi);

// Check the moduli interpretation
jj,jj_1728 := WNrModuli(14,1 : vars:=[(1 + e^2*t)/e, 1/e]);
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
table := eval Read("../Z14-r/curvesonZ14-r/modularcurves/14-1.m");

for ex in table do
  m, X := Explode(ex);
  printf "The case when m=%o\n", m;
  X := Curve(ZZ, [X]);
  X0m := SmallModularCurve(m);
  g := Genus(X0m);
  assert Genus(X) eq g;

  if g eq 0 then
    KK := FunctionField(QQ);
    coords := [KK.1];
    gen_pt := [KK.1, 1];
    X0m := BaseChange(X0m, KK);
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

  jj, jj_1728 := WNrModuli(14,1 : vars:=pp);
  j1 := jNInvariant(X0m!gen_pt, m);
  j2 := jInvariant(X0m!gen_pt, m);

  assert jj eq j1*j2;
  assert jj_1728 eq (j1 - 1728)*(j2 - 1728);
end for;


/**************************************************
  Z(14,3)
*/
ZZ := Surface(A3, ZNrEquations(14,3 : vars:=[u,v,z]));
print "--------------------\nThe case Z(14,3)\n";

////////////////////
// FIRST DEAL WITH THE BLOWN DOWN CASES
print "The blown down cases:\n";
f := -ZNrEquations(14,3 :  vars:=[x[1]/x[3],x[2]/x[3],0])[1];
f := x[3]^14 * f;

// X_0(5)
printf "The case when m=5\n";
bl_up := Evaluate(f, [1/2*e, 1 + 1/4*e^2*t, 1]);
assert bl_up + BigO(e^5) eq (-16*t^2 + 88*t + 4)*e^4  + BigO(e^5);

// Parametrise the exceptional divisor (where e=0)
pp := [
    -2*t/(t^2 + 22*t + 125),
    (-2*t^2 + 250)/(t^2 + 22*t + 125)
];

F_5 := Curve(A2, b^2 - (-16*a^2 + 88*a + 4));
psi := map<A1 -> F_5 | [Evaluate(p, A1.1) : p in pp]>;
assert IsInvertible(psi);

// Check the moduli interpretation
jj,jj_1728 := WNrModuli(14,3 : vars:=[1/2*e, 1 + 1/4*e^2*t]);
jj := Evaluate(jj,0); jj_1728 := Evaluate(jj_1728,0); // On the exceptional divisor
jj := Evaluate(jj, pp[1]); jj_1728 := Evaluate(jj_1728, pp[1]);

X0m := SmallModularCurve(5);
X0m := BaseChange(X0m, Kt);
j1 := jNInvariant(X0m![t,1], 5);
j2 := jInvariant(X0m![t,1], 5);

assert j1*j2 eq jj;
assert (j1-1728)*(j2-1728) eq jj_1728;

//////////////////////////////////////////////////.
// X_0(13)
printf "The case when m=13\n";
bl_up := Evaluate(f, [-4 + 1/4*e, -1 + 1/12*e + 1/4*e*t, 1]);
assert bl_up + BigO(e^5) eq -t^2*(477*t^2 + 48*t - 1)*e^4  + BigO(e^5);

// Parametrise the exceptional divisor (where e=0)
pp := [
    2/9*t/(t^2 + 16/3*t + 13),
    (t^2 - 13)/(t^2 + 16/3*t + 13)
];

F_13 := Curve(A2, b^2 + (477*a^2 + 48*a - 1));
psi := map<A1 -> F_13 | [Evaluate(p, A1.1) : p in pp]>;
assert IsInvertible(psi);

// Check the moduli interpretation
jj,jj_1728 := WNrModuli(14,3 : vars:=[-4 + 1/4*e, -1 + 1/12*e + 1/4*e*t]);
jj := Evaluate(jj,0); jj_1728 := Evaluate(jj_1728,0); // On the exceptional divisor
jj := Evaluate(jj, pp[1]); jj_1728 := Evaluate(jj_1728, pp[1]);

X0m := SmallModularCurve(13);
X0m := BaseChange(X0m, Kt);
j1 := jNInvariant(X0m![t,1], 13);
j2 := jInvariant(X0m![t,1], 13);

assert j1*j2 eq jj;
assert (j1-1728)*(j2-1728) eq jj_1728;

//////////////////////////////////////////////////.
// X_0(17)
printf "The case when m=17\n";
bl_up := Evaluate(f, [1/4*e, -1 + 1/4*e*t, 1]);
assert bl_up + BigO(e^5) eq -(t^4 - 16*t^3 + 9*t^2 - 2*t + 4)*e^4  + BigO(e^5);

jj,jj_1728 := WNrModuli(14,3 : vars:=[1/4*e, -1 + 1/4*e*t]);
jj := Evaluate(jj,0); jj_1728 := Evaluate(jj_1728,0); // On the exceptional divisor

X0m := SmallModularCurve(17);
KK<x0,x1> := FunctionField(X0m);

// Check that this curve is isomorphic to the exceptional divisor
pp := [
  -2/(x0^2 - 6*x0 + 10)*x1 + (x0^2 - 6)/(x0^2 - 6*x0 + 10),
  (14*x0^2 - 60*x0 + 40)/(x0^4 - 12*x0^3 + 56*x0^2 - 120*x0 + 100)*x1 + (-2*x0^4 - 2*x0^3 + 40*x0^2 + 64*x0 - 320)/(x0^4 - 12*x0^3 + 56*x0^2 - 120*x0 + 100)
];
F_17 := Curve(A2, b^2 + (a^4 - 16*a^3 + 9*a^2 - 2*a + 4));
psi := map<X0m -> F_17 | pp>;
assert IsInvertible(psi);

// Check the moduli interpretation
jj := Evaluate(jj, pp[1]); jj_1728 := Evaluate(jj_1728, pp[1]);

X0m_KK := BaseChange(X0m, KK);
j1 := jNInvariant(X0m_KK![x0,x1,1], 17);
j2 := jInvariant(X0m_KK![x0,x1,1], 17);

assert j1*j2 eq jj;
assert (j1-1728)*(j2-1728) eq jj_1728;


////////////////////
// NOW THE NON-BLOWN DOWN
print "\nAnd the remaining cases:\n";
table := eval Read("../Z14-r/curvesonZ14-r/modularcurves/14-3.m");

for ex in table do
  m, X := Explode(ex);
  printf "The case when m=%o\n", m;
  X := Curve(ZZ, [X]);
  X0m := SmallModularCurve(m);
  g := Genus(X0m);
  assert Genus(X) eq g;

  if g eq 0 then
    KK := FunctionField(QQ);
    coords := [KK.1];
    gen_pt := [KK.1, 1];
    X0m := BaseChange(X0m, KK);
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

  jj, jj_1728 := WNrModuli(14,3 : vars:=pp);
  j1 := jNInvariant(X0m!gen_pt, m);
  j2 := jInvariant(X0m!gen_pt, m);

  assert jj eq j1*j2;
  assert jj_1728 eq (j1 - 1728)*(j2 - 1728);
end for;

