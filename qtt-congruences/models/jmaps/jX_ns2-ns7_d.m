X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns2-ns7_d.m");

x := X2.1/X2.3; y := X2.2/X2.3;
t := -x;
phi := map<X2 -> X1 | [((3*t + 1)*(t^2 + 3*t + 4)*(t^2 + 10*t + 4)*(4*t^2 + 5*t + 2))^3/(t^3 + t^2 - 2*t - 1)^7, 1]>;
return phi;