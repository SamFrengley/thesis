P2<X,Y,Z> := ProjectiveSpace(Rationals(), 2);
fs := [X^3*Y + X^2*Y^2 + X*Y^3 + 3*X^3*Z + 2*X*Y^2*Z + Y^3*Z + 3*X*Y*Z^2 + Y^2*Z^2 
    + Y*Z^3];
return Curve(P2, fs);

