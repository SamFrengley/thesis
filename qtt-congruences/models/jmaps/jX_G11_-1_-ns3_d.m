X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_G11_-1_-ns3_d.m");

x := X2.1/X2.3; y := X2.2/X2.3;
s := (-2*(3*x + y + 15)/(3*x - y + 15));
t3 := (s^6 - 16)/(s^2);
phi := map<X2 -> X1 | [t3^3, 1]>;
return phi;