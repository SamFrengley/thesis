X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns4-ns3_d.m");

t := X2.1/X2.2;
t4 := (t^2 - 2*t - 2)^3/(24*(t^3 -4)*(t^3+2);
phi := map<X2 -> X1 | [(32*t4 - 4)/t4^4]>;
return phi;