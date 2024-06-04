X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns2-ns5_d.m");

t := X2.1/X2.2;
phi := map<X2 -> X1 | [(5^3*(t + 1)*(2*t + 1)^3*(2*t^2 - 3*t + 3)^3)/((t^2 + t - 1)^5), 1]>;
return phi;