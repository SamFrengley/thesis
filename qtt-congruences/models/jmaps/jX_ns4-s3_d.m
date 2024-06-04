X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns4-s3_d.m");

x := X2.1/X2.3; y := X2.2/X2.3;
s := (6*x^3 + x^2*y + 14*x^2 + 12*x*y - 8*x + 8*y - 16)/((x + 4)^2*(2*x + y + 2));
t3 := 3*(s^2-1)*(3*s^2+1)/(s^2);
phi := map<X2 -> X1 | [t3^3, 1]>;
return phi;