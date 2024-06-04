P3<X,Y,Z,W> := ProjectiveSpace(Rationals(), 3);
fs := [
    Y^2*Z - Y*Z^2 + 2*Y^2*W + 4*Y*Z*W - 4*Z^2*W + 12*Y*W^2 + 23*W^3,
    X^2 + X*Y + Y*Z - Z^2 - X*W + Y*W + W^2
];
return Curve(P3, fs);
