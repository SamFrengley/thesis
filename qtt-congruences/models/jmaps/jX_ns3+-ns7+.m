X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns3+-ns7+.m");

x := X2.1/X2.3; y := X2.2/X2.3;

t := (6 *x^2 - x *y - 11 *x  - 7 *y  + 70 )/(3*x^2 - x *y + 10 *x  - 7 *y  - 77);
j := ((3*t + 1)*(t^2 + 3*t + 4)*(t^2 + 10*t + 4)*(4*t^2 + 5*t + 2))^3/(t^3 + t^2 - 2*t - 1)^7;
phi := map<X2 -> X1 | [j, 1]>;
return phi;