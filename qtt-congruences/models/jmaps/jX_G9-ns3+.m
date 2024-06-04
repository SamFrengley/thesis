X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_G9-ns3+.m");

x := X2.1/X2.3; y := X2.2/X2.3;

phi := map<X2 -> X1 | [(4096*y^6 + 9216*y^4 + 6912*y^2 + 1728)/(y^2 + 1), 1]>;
return phi;