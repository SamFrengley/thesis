X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns2-s5_d.m");

x := X2.1/X2.3; y := X2.2/X2.3;
phi := map<X2 -> X1 | [((x + 3)*(x^2 + x + 4)*(x^2- 4*x - 1))^3/(x^2 + x - 1)^5, 1]>;
return phi;