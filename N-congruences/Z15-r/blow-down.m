
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//(15,2)

Attach("../ZNr-equations.m");
QQ := Rationals();
A3<u,v,z> := AffineSpace(QQ, 3);
_<uu,vv> := FunctionField(QQ,2);
jj,jj_1728 := WNrModuli(15,11 : vars:=[1/uu,vv/uu]);

//f := ZNrEquations(15,2:vars:=[u,v,z])[1];
f := -ZNrEquations(15,11:vars:=[1/u,v/u,0])[1];
f := z^2 - Numerator(f);

S := Surface(A3, f); SingularPoints(S);
//pts := SingularSS(S); pts;
//pts := [Scheme(S, [u,v,z]), Scheme(S, [u+1,v,z])];

pts := [Scheme(S, [u, v,z])];
dsd := ResolveSingByBlowUp(S, pts[1]);

n :=  NumberOfBlowUpDivisors(dsd);
printf "The number of blowup divisors is %o\n\n", n;

PrintFile("store2.m", "[");
for i in [1..n] do
  K<e,t,z> := FunctionField(QQ,3);
  E,_,m := BlowUpDivisor(S,dsd,i);
  pp := DefiningPolynomials(m)[1..2];
  E := DefiningPolynomials(E);
  pp := [Evaluate(p, [e,t,z]) : p in pp];
  E := [Evaluate(p, [e,t,z]) : p in E];
  PrintFile("store2.m", <E,pp>);
  if i ne n then
    PrintFile("store2.m", ",");
  end if;
end for;

PrintFile("store2.m", "]");


/* n :=  NumberOfBlowUpDivisors(dsd); */
/* printf "The number of blowup divisors is %o\n\n", n; */

/* for i in [1..n] do */
/*   try */
/*     E,_,m := BlowUpDivisor(S, dsd, i); */
/*     pp := DefiningPolynomials(m)[1..2]; */

/*     printf "The genus of the exceptional is %o\n\n", Genus(E); */

/*     j := Evaluate(jj,pp); */
/*     j_1 := Evaluate(jj_1728,pp); */
/*     j_s := (j - j_1 + 1728^2)/1728; */


/*     pp := PointSearch(E, 80); */
/*     printf "We found %o points", #pp; */
/*     p := pp[Random([1..#pp])]; */
/*     j := Evaluate(j, Eltseq(p)); */
/*     js := Evaluate(j_s, Eltseq(p)); */


/*     _<x> := PolynomialRing(QQ); */
/*     g := x^2 - js*x + j; */
/*     j1,j2 := Explode([r[1] : r in Roots(g)]); */
/*     E1 := EllipticCurveWithjInvariant(j1); */
/*     IsogenousCurves(E1); */
/*   catch e */
/*     e; */
/*   end try; */
/* end for; */
  
//////////////////////////////////////////////////.
Kt<t> := FunctionField(QQ);
P<e> := PowerSeriesRing(Kt);
A2<x,y> := AffineSpace(Rationals(), 2);
_<X> := PolynomialRing(Kt);


