X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns2-ns9_d.m");

x := X2.1/X2.3; y := X2.2/X2.3;
t := (3-x)/x;
t3 := -3*(t^3 + 3*t^2 - 6*t + 4)*(t^3 + 3*t^2 + 3*t + 4)*(5*t^3 - 3*t^2 - 3*t + 2)/(t^3 - 3*t + 1)^3;
phi := map<X2 -> X1 | [t3^3, 1]>;
return phi;