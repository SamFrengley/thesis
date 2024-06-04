Attach("../ZNr-equations.m");
/* QQ := Rationals(); */
/* A3<u,v,z> := AffineSpace(QQ, 3); */
/* jj,jj_1728 := WNrModuli(14,1); */
/* P := Parent(jj); */
/* jj := Evaluate(jj, [P.1/P.2, 1/P.2]); */
/* jj_1728 := Evaluate(jj_1728, [P.1/P.2, 1/P.2]); */

/* f := ZNrEquations(14,1)[1]; */
/* f := -Evaluate(f, [u/v,1/v,0]); */
/* f := Numerator(f); */
/* f := z^2 - f; */

/* S := Surface(A3, f); */
/* pts := SingularSS(S); */
/* dsd := ResolveSingByBlowUp(S, pts[2]); */
/* printf "The number of blowup divisors is %o\n\n", NumberOfBlowUpDivisors(dsd); */
/* E,_,m := BlowUpDivisor(S, dsd, 3); */
/* pp := DefiningPolynomials(m)[1..2]; */

/* printf "The genus of the exceptional is %o\n\n", Genus(E); */

/* jj := Evaluate(jj,pp); */
/* jj_1728 := Evaluate(jj_1728,pp); */
/* j_sum := (jj - jj_1728 + 1728^2)/1728; */


/* pp := PointSearch(E, 100); #pp; */
/* p := pp[Floor(#pp/2)]; */
/* j := Evaluate(jj, Eltseq(p)); */
/* js := Evaluate(j_sum, Eltseq(p)); */


/* _<x> := PolynomialRing(QQ); */
/* f := x^2 - js*x + j; */
/* j1,j2 := Explode([r[1] : r in Roots(f)]); */
/* E1 := EllipticCurveWithjInvariant(j1); */
/* IsogenousCurves(E1); */

/* //////////////////////////////////////////////////. */
/* Kt<t> := FunctionField(QQ); */
/* P<e> := PowerSeriesRing(Kt); */
/* A2<x,y> := AffineSpace(Rationals(), 2); */
/* _<X> := PolynomialRing(Kt); */

/* /\* f := ZNrEquations(14,1)[1]; *\/ */
/* /\* jj,jj_1728 := WNrModuli(14,1); *\/ */
/* /\* j_sum := (jj - jj_1728 + 1728^2)/1728; *\/ */

/* /\* ////////////////////////////////////////////////// *\/ */
/* /\* // X_0(9) *\/ */
/* /\* bl_up := Evaluate(-f, [e,e^2*t - e,0]); *\/ */
/* /\* assert bl_up + BigO(e^5) eq (t^2 - 16*t - 44)*e^4 + BigO(e^5); *\/ */
/* /\* E := Curve(A2, y^2 - (x^2 - 16*x - 44)); *\/ */
/* /\* E := ProjectiveClosure(E); *\/ */
/* /\* /\\* pp := DefiningPolynomials(Parametrization(E, E![1,1,0])); *\\/ *\/ */
/* /\* /\\* pp := [Evaluate(p, [t,1]) : p in pp]; *\\/ *\/ */
/* /\* /\\* pp := [pp[1]/pp[3], pp[2]/pp[3]]; *\\/ *\/ */
/* /\* /\\* pp := [Evaluate(p, -2*t - 8) : p in pp]; *\\/ *\/ */
/* /\* pp := [ *\/ */
/* /\*   (t^2 + 8*t + 27)/t, *\/ */
/* /\*   (t^2 - 27)/t *\/ */
/* /\* ]; *\/ */

/* /\* jj := Evaluate(jj, [e,e^2*t - e]); *\/ */
/* /\* jj_1728 := Evaluate(jj_1728, [e,e^2*t - e]); *\/ */
/* /\* jj := Evaluate(jj,0); *\/ */
/* /\* jj_1728 := Evaluate(jj_1728,0); *\/ */
/* /\* jj := Evaluate(jj, pp[1]); *\/ */
/* /\* jj_1728 := Evaluate(jj_1728, pp[1]); *\/ */

/* /\* X0m := SmallModularCurve(9); *\/ */
/* /\* X0m := BaseChange(X0m, Kt); *\/ */
/* /\* j1 := jNInvariant(X0m![t,1], 9); *\/ */
/* /\* j2 := jInvariant(X0m![t,1], 9); *\/ */
/* /\* assert j1*j2 eq jj; *\/ */
/* /\* assert (j1-1728)*(j2-1728) eq jj_1728; *\/ */

/* /\* /\\* js := Evaluate(j_sum, [e,e^2*t - e]); *\\/ *\/ */
/* /\* /\\* j := Evaluate(j, 0); *\\/ *\/ */
/* /\* /\\* js := Evaluate(js, 0); *\\/ *\/ */
/* /\* /\\* j := Evaluate(j, pp[1]); *\\/ *\/ */
/* /\* /\\* js := Evaluate(js, pp[1]); *\\/ *\/ */
/* /\* /\\* f := X^2 - js*X + j; *\\/ *\/ */
/* /\* /\\* j1,j2 := Explode([r[1] : r in Roots(f)]); *\\/ *\/ */

/* /\* /\\* X0m := SmallModularCurve(9); *\\/ *\/ */
/* /\* /\\* X0m := BaseChange(X0m, Kt); *\\/ *\/ */
/* /\* /\\* j1 := jNInvariant(X0m![t,1], 9); *\\/ *\/ */
/* /\* /\\* j2 := jInvariant(X0m![t,1], 9); *\\/ *\/ */

/* /\* /\\* Roots(Numerator(Evaluate(j1*j2, t) - Evaluate(j, X))); *\\/ *\/ */



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//(14,3)

Attach("../ZNr-equations.m");
QQ := Rationals();
A3<u,v,z> := AffineSpace(QQ, 3);
_<uu,vv> := FunctionField(QQ,2);
jj,jj_1728 := WNrModuli(14,1 : vars:=[1/uu,vv/uu]);
u := 1/u; v := v/u;

f := -ZNrEquations(14,1:vars:=[u,v,0])[1];

f := z^2 - Numerator(f);
S := Surface(A3, f);

/* S2 := Scheme(S, [Numerator(v-1), Numerator(u*v^2 + 4*u*v - u - v^3 - 5*v^2 + 13*v + 1)]); */
/* pts := IrreducibleComponents(ReducedSubscheme(S2)); */
pts := SingularSS(S);
pts;
dsd := ResolveSingByBlowUp(S, pts[5]); 

n :=  NumberOfBlowUpDivisors(dsd);
printf "The number of blowup divisors is %o\n\n", n;

/* PrintFile("store.m", "["); */
/* for i in [1..n] do */
/*   K<e,t,z> := FunctionField(QQ,3); */
/*   E,_,m := BlowUpDivisor(S,dsd,i); */
/*   pp := DefiningPolynomials(m)[1..2]; */
/*   E := DefiningPolynomials(E); */
/*   pp := [Evaluate(p, [e,t,z]) : p in pp]; */
/*   E := [Evaluate(p, [e,t,z]) : p in E]; */
/*   PrintFile("store.m", <E,pp>); */
/*   if i ne n then  */
/*     PrintFile("store.m", ","); */
/*   end if; */
/* end for; */

/* PrintFile("store.m", "]"); */

for i in [1..n] do
  try
    E,_,m := BlowUpDivisor(S, dsd, i);
    pp := DefiningPolynomials(m)[1..2];

    printf "The genus of the exceptional is %o\n\n", Genus(E);

    j := Evaluate(jj,pp);
    j_1 := Evaluate(jj_1728,pp);
    j_s := (j - j_1 + 1728^2)/1728;


    pp := PointSearch(E, 40);
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
  
/* //////////////////////////////////////////////////. */
/* Kt<t> := FunctionField(QQ); */
/* P<e> := PowerSeriesRing(Kt); */
/* A2<x,y> := AffineSpace(Rationals(), 2); */
/* _<X> := PolynomialRing(Kt); */

/* A3<e,t,z> := AffineSpace(QQ,3); */
/* bud := eval Read("store.m"); n:= #bud; */


/* function FactEval(f, pp) */
/*   f_n := Factorisation(Numerator(f)); */
/*   f_d := Factorisation(Denominator(f)); */
/*   f2 := (&*[ff[1]^ff[2] : ff in f_n])/(&*[ff[1]^ff[2] : ff in f_d]); */
/*   c := Rationals()!(f2/f); */
  
/*   f_n := [<Evaluate(ff[1], pp), ff[2]> : ff in f_n]; */
/*   f_d := [<Evaluate(ff[1], pp), ff[2]> : ff in f_d]; */
/*   f2 := (&*[ff[1]^ff[2] : ff in f_n])/(&*[ff[1]^ff[2] : ff in f_d]); */
/*   return c*f2; */
/* end function; */
  
/* Kt<t> := FunctionField(QQ); */
/* P<e> := PowerSeriesRing(Kt); */


/* bud := Reverse(bud); */
/* for i in [1..n] do */
/*   try */
/*     E := Curve(A3, bud[i][1]); */
/*     E1 := ProjectiveClosure(E); */
/*     phi1 := map<E -> E1 | [E.1,E.2,E.3,1]>; */
/*     phi2 := map<E1 -> ProjectiveSpace(QQ,2) | Basis(-CanonicalDivisor(E1))>; */
/*     E2 := Image(phi2); phi2 := map<E1 -> E2 | DefiningPolynomials(phi2)>; */
/*     _, PHI := IsInvertible(phi1*phi2); */
/*     points_E2 := PointSearch(E2, 50000); */
/*     points := []; */
/*     for p in points_E2 do */
/*       try */
/*         Append(~points, PHI(p)); */
/*       catch e */
/*         assert 0 eq 0; */
/*       end try; */
/*     end for; */

/*     printf "The genus of the exceptional is %o\n", Genus(E); */
/*     printf "We found %o points\n\n", #points; */

/*     pp := bud[i][2]; */

/*     j := FactEval(jj, pp); */
/*     j_1 := FactEval(jj_1728, pp); */
/*     j_s := (j - j_1 + 1728^2)/1728; */

/*     p := points[Random([1..#points])]; */
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
