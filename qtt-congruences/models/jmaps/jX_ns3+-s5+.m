X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns3+-s5+.m");

x := X2.1/X2.3; y := X2.2/X2.3;

t3 := - x *(y - 2)*(y^2 + y + 4)*(y^2 + 6*y + 4)/(y^2 + y - 1)^2;
phi := map<X2 -> X1 | [t3^3, 1]>;
return phi;