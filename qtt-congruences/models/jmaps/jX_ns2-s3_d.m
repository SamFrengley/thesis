X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns2-s3_d.m");

s := X2.1/X2.2;
t3 := 3*(s^2-1)*(3*s^2+1)/(s^2)
phi := map<X2 -> X1 | [t3^3, 1]>;
return phi;