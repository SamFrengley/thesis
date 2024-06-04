/*
    Functions for determining all the 4 and 8 congruent quadratic 
    twists of an elliptic curve. The method is via the equations for X_E^r(N).

    Equations and forgetful maps for X_E^r(8) are taken from Tom Fisher's code, see 
    https://www.dpmms.cam.ac.uk/~taf1000/papers/congr-ellsurf.html

    Functions include 

    X8ErEquations(ab, r)
        // If E : y^2 = x^3 + ax + b, returns the curve X_E^r(8) 
        // and the map X_E^r(8) --> X_E^r(4) given by Chen.

    NonTrivially4CongruentQT(ab, r)
        // Let E: y^2 = x^3 + ax + b
        // Returns the quadratic twists of E which are nontrivially (4, r)-congruent 
        // to E and the corrseponding points on Fisher's X_E^r(4).

    NonTrivially8CongruentQT(ab, r)
        // Let E: y^2 = x^3 + ax + b
        // Returns the quadratic twists of E which are nontrivially (8, r)-congruent 
        // to E and the corrseponding points on Chen's X_E^r(8).
*/

// Equations for X_E^r(8) and forgetful maps

function X8ErEquations(ab, r)
    // If E : y^2 = x^3 + ax + b, returns the curve X_E^r(8) 
    // and the map X_E^r(8) --> X_E^r(4) given by Chen.
    P1 := ProjectiveSpace(Rationals(), 1);
    a,b := Explode(ab);
    case r : 
    when 1 :
        PP<x0,x1,x2,x3,x4> := ProjectiveSpace(Rationals(),4);
        // Equations for X_E(8,1) --- Chen
        FF := [ -a*x3^2 + 2*x1*x3 + x2^2 + 2*x4^2,
          -2*a*x2*x3 - b*x3^2 + 2*x1*x2 + 2*x0*x4,
          -2*b*x2*x3 + x1^2 - x0^2 + a*x4^2];
        XE := Curve(PP, FF);
        s := x0;
        t := 3*x4;
        phi := map<XE -> P1 | [s,t]>;
    when 3 :
        PP<x0,x1,x2,x3,x4> := ProjectiveSpace(Rationals(),4);
        // Equations for X_E(8,3) --- Chen
        FF := [ x0^2 + 2*b*x1*x3 + a*x2^2 - 6*b*x2*x4 - a^2*x4^2,
           2*x0*x1 + 2*a*x1*x3 + 4*a*x2*x4 - b*x3^2 - 18*b*x4^2,
           2*x0*x3 - x1^2 + x2^2 + a*x3^2 + 3*a*x4^2 ];
        XE := Curve(PP, FF);
        s := (-3*b*x2 - 2*a^2*x4);
        t := (6*a*x2 - 27*b*x4);
        phi := map<XE -> P1 | [s,t]>;
    when 5 : 
        PP<x0,x1,x2,x3,x4> := ProjectiveSpace(Rationals(),4);
        // Equations for X_E(8,5) --- Chen
        FF := [ -a*x3^2 + 2*x1*x3 + x2^2 - (5*a*x4^2 - 3*x0^2),
           -2*a*x2*x3 - b*x3^2 + 2*x1*x2 - (2*a*x0*x4 - 6*b*x4^2),
           -2*b*x2*x3 + x1^2 - (4*a^2*x4^2 - 4*a*x0^2 - 6*b*x0*x4) ];
        XE := Curve(PP, FF);
        s := x0;
        t := 3*x4;
        phi := map<XE -> P1 | [s,t]>;
    when 7 :
        PP<x0,x1,x2,x3,x4> := ProjectiveSpace(Rationals(),4);
        // Equations for X_E(8,7) --- Chen
        FF := [ 3*x0^2 + a*x4^2 - a*x3^2 + 2*x1*x3 + x2^2,
           4*a*x0*x4 + 6*b*x4^2 - 2*a*x2*x3 - b*x3^2 + 2*x1*x2,
           a*x0^2 + 6*b*x0*x4 - a^2*x4^2 - 2*b*x2*x3 + x1^2 ];
        XE := Curve(PP, FF);
        s := x0;
        t := 3*x4;
        phi := map<XE -> P1 | [s,t]>;
    end case;
    return XE, phi;
end function;

// Functions for determining congruent quadratic twists.

function NonTrivially4CongruentQT(ab, r)
    // Let E: y^2 = x^3 + ax + b
    // Returns the quadratic twists of E which are nontrivially (4, r)-congruent 
    // to E and the corrseponding points on Fisher's X_E^r(4).

    XEr := ProjectiveSpace(Rationals(), 1);
    E := EllipticCurve(ab);
    cc := [-ab[1]/27, -ab[2]/54];
    j := jInvariant(E);
    Delta := Discriminant(E);

    _, c4, c6 := HessePolynomials(4, 1, cc : Variables:=[XEr.1,XEr.2]);
    f := (1728 - j)*c4^3 + j*c6^2;
    pts := Points(Cluster(XEr, f));

    ret := [];

    for pp in pts do
        p := Coordinates(pp);
        if r eq 1 then
            Ed := EllipticCurve([-27*Evaluate(c4, p), -54*Evaluate(c6, p)]);
        elif r eq 3 then 
            Ed := EllipticCurve([-27*Delta^2*Evaluate(c4, p), -54*Delta^3*Evaluate(c6, p)]);
        end if;
        Ed := MinimalModel(Ed);

        if IsQuadraticTwist(E, Ed) then
            flag, D := IsIsogenous(E,Ed);
            if (not flag) then
                Append(~ret, <Ed, p>);
            end if;
        end if;
    end for;

    return ret;
end function;

function NonTrivially8CongruentQT(ab, r)
    // Let E: y^2 = x^3 + ax + b
    // Returns the quadratic twists of E which are nontrivially (8, r)-congruent 
    // to E and the corrseponding points on Chen's X_E^r(8).
    
    XE8, phi := X8ErEquations(ab, r);
    XE4 := Codomain(phi);

    pts4 := NonTrivially4CongruentQT(ab, r mod 4);
    
    ret := [];
    for p in pts4 do
        pts8 := Points((XE4!p[2]) @@ phi);
        if #pts8 gt 0 then
            ret := ret cat [<p[1], Coordinates(pp)> : pp in pts8];
        end if;
    end for;

    return ret;

end function;

function TrFrobDiff(E1, E2, p)
    return TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p);
end function;
