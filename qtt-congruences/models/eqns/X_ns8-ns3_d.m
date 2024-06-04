P2<X,Y,Z> := ProjectiveSpace(Rationals(), 2);
fs := [X^3*Y + 1/2*X^2*Y^2 - X*Y^3 - 1/2*Y^4 - X^3*Z - X^2*Y*Z + X*Y^2*Z + Y^3*Z - 
   X^2*Z^2 + X*Y*Z^2 + 1/2*Y^2*Z^2 - X*Z^3 - Y*Z^3];
return Curve(P2, fs);