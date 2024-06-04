P<X,Y,Z,W> := ProjectiveSpace(Rationals(), 3);
fs := [
    X*Y^2 + X^2*Z - 3*X*Y*Z - Y*Z^2 + X*W^2 - 2*Y*W^2,
    X^2 + 2*X*Y - 2*X*Z + Y*Z - Z^2 - W^2
];
return Curve(P, fs);