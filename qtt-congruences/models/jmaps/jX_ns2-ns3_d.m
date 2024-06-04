X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns2-ns3_d.m");

t := X2.1/X2.2;
phi := map<X2 -> X1 | [1728*(1-t^2)^3), 1]>;
return phi;