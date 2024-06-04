
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//(15,2)

Attach("../ZNr-equations.m");
QQ := Rationals();
A3<u,v,z> := AffineSpace(QQ, 3);
_<uu,vv> := FunctionField(QQ,2);
jj,jj_1728 := WNrModuli(15,2 : vars:=[uu,vv]);

h := Numerator(jj);
k := Denominator(jj);
h := &*[a[1] : a in FactNum(h)];
k := &*[a[1] : a in FactNum(k)];

//f := ZNrEquations(15,2:vars:=[u,v,z])[1];
f := z^2 + Numerator(ZNrEquations(15,2:vars:=[u,v,0])[1]);
S := Surface(A3, f);

pts := [Scheme(S, [u,v,z]), Scheme(S, [u,v-1/2,z]), Scheme(S, [u-2,v-2,z])];
pt := pts[2];
cmps := LocalBlowUp(S, pt); n := #cmps;
printf "The number of blowup divisors is %o\n\n", n;

for i in [1..n] do
  try
    m := cmps[i][2];
    pp := DefiningPolynomials(m)[1..2];
    E := IrreducibleComponents((pt @@ m));
    E := Curve(E[1]);   
      
    printf "The genus of the exceptional is %o\n\n", Genus(E);

    j := Evaluate(jj,pp);
    j_1 := Evaluate(jj_1728,pp);
    j_s := (j - j_1 + 1728^2)/1728;


    pp := PointSearch(E, 80);
    printf "We found %o points", #pp;
    p := pp[Random([1..#pp])];
    j := Evaluate(j, Eltseq(p));
    js := Evaluate(j_s, Eltseq(p));


    _<x> := PolynomialRing(QQ);
    g := x^2 - js*x + j;
    j1,j2 := Explode([r[1] : r in Roots(g)]);
    E1 := EllipticCurveWithjInvariant(j1);
    IsogenousCurves(E1);
  catch e
    e;
  end try;
end for;
  
//////////////////////////////////////////////////.
Kt<t> := FunctionField(QQ);
P<e> := PowerSeriesRing(Kt);
A2<x,y> := AffineSpace(Rationals(), 2);
_<X> := PolynomialRing(Kt);


