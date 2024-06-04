/*
    We claim that some maps between curves are birational in Sections 4.4 and 4.5.
    The code in this file checks this. 
*/
QQ := Rationals();

// ----------------
// Lemma 4.3
// ----------------
A2<s,u> := AffineSpace(QQ, 2);
f := 4*s^2*u*(u^3 - 4) - (-3*(s^2 - 1)*(3*s^2 + 1));
C1 := Curve(A2, f);

C2 := EllipticCurve([0,1,0,-4,-4]);
_<x,y> := FunctionField(C2);

su := [(6*x^3 + x^2* y + 14* x^2 + 12* x* y - 8* x + 8* y - 16)/((x + 4)^2 *(2 *x + y + 2)) , 6*(x^2 + 3 *x + y + 2)/ ((x + 4) *(2 *x + y + 2))];
assert IsInvertible(map<C2 -> C1 | su>);

// ----------------
// Lemma 4.8
// ----------------
A2<s,u> := AffineSpace(QQ, 2);
f := u^2 - 3*(s^2 - 2*s + 4)*(s^2 + 2*s + 4);
C1 := Curve(A2, f);

C2 := EllipticCurve([-39,-70]);
_<x,y> := FunctionField(C2);

su := [-2*(3*x + y + 15)/(3*x - y + 15), 12*(x - 1)*(x + 5)*(x + 11)/(3*x - y + 15)^2];
assert IsInvertible(map<C2 -> C1 | su>);

// ----------------
// Lemma 5.7
// ----------------
A3<u,v,w> := AffineSpace(QQ, 3);
C1 := Curve(A3, [v^2 + v - u^3 - 1, w^2 - 3*(3*u^2 -2*u - 4*v + 3)*(9*u^2 - 10*u + 5)]);

C2 := eval Read("../models/eqns/X_ns3-ns5_d.m"); 
_<x,y> := FunctionField(C2);

uvw := [(x^4 + 2*x^3 - x + 3*y + 4)/(x^4 + 2*x^3 - 3*x^2 - 4*x + 4),
        (x^6 + 3*x^5 + 6*x^4 + 7*x^3 + 3*x^2*y - 6*x^2 + 3*x*y - 9*x + 3*y + 
            16)/(x^6 + 3*x^5 - 3*x^4 - 11*x^3 + 6*x^2 + 12*x - 8),
        (6*x^7 + 21*x^6 + 135*x^5 + 285*x^4 + 36*x^3*y - 51*x^3 + 54*x^2*y - 351*x^2
            + 108*x*y + 249*x + 45*y + 192)/(x^8 + 4*x^7 - 2*x^6 - 20*x^5 + x^4 + 
            40*x^3 - 8*x^2 - 32*x + 16)
];
assert IsInvertible(map<C2 -> C1 | uvw>);

// ----------------
// Lemma 5.12
// ----------------
A3<u,v,w> := AffineSpace(QQ, 3);
C1 := Curve(A3, [v^2 + v - u^3 + u^2 + 7*u - 10, w^2 - ((-3*u^2 + 24*u - 37)*v - u^4 + 5*u^3 + 18*u^2 - 95*u + 94)]);
C2<X,Y,Z,W> := eval Read("../models/eqns/X_ns2-ns11_d.m");

uvw := [(-18/11*X^2*Y*Z + 6/11*X^2*Z^2 - 391/55*X^2*W^2 - 28/11*X*Y^2*Z + 
            102/11*X*Y*Z^2 - 299/55*X*Y*W^2 - 48/11*X*Z^3 + 102/11*X*Z*W^2 + 
            46/11*Y^3*Z - 108/11*Y^2*Z^2 + 138/11*Y^2*W^2 + 48/11*Y*Z^3 - 
            458/55*Y*Z*W^2 + 8/11*Z^4 + 12/5*Z^2*W^2 + 92/55*W^4)/(X*Z*W^2 - 
            42/55*Y*Z*W^2 + 2/11*Z^4 + 3/5*Z^2*W^2 + 23/55*W^4),
        (36/11*X^2*Y*Z - 12/11*X^2*Z^2 + 782/55*X^2*W^2 + 38/11*X*Y^2*Z - 18*X*Y*Z^2
            + 207/55*X*Y*W^2 + 96/11*X*Z^3 - 182/11*X*Z*W^2 - 138/11*Y^3*Z + 
            324/11*Y^2*Z^2 - 414/11*Y^2*W^2 - 144/11*Y*Z^3 + 102/5*Y*Z*W^2 - 
            12/11*Z^4 - 18/5*Z^2*W^2 - 138/55*W^4)/(X*Z*W^2 - 42/55*Y*Z*W^2 + 
            2/11*Z^4 + 3/5*Z^2*W^2 + 23/55*W^4),
        (18/11*X*Y*Z*W - 6/11*X*Z^2*W + 391/55*X*W^3 + 46/11*Y^2*Z*W - 
            108/11*Y*Z^2*W + 138/11*Y*W^3 + 48/11*Z^3*W - 58/11*Z*W^3)/(X*Z*W^2 - 
            42/55*Y*Z*W^2 + 2/11*Z^4 + 3/5*Z^2*W^2 + 23/55*W^4)
];
assert IsInvertible(map<C2 -> C1 | uvw>);
