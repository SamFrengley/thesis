X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_s4+-ns3+.m");

x := X2.1/X2.3; y := X2.2/X2.3;
j := (-64*y^12 + 384*y^10 - 192*y^8 - 1792*y^6 + 576*y^4 + 3456*y^2 + 1728)/(y^8 - 4*y^6 + 6*y^4 - 4*y^2 + 1);
phi := map<X2 -> X1 | [j, 1]>;
return phi;