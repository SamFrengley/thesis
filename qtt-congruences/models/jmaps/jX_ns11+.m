X1 := eval Read("../models/eqns/X1.m");
X2 := eval Read("../models/eqns/X_ns11+.m");

x := X2.1/X2.3; y := X2.2/X2.3;

f1 := x^2 + 3*x - 6;
f2 := 11*(x^2 - 5)*y + (2*x^4 + 23*x^3 - 72*x^2 - 28*x + 127);
f3 := 6*y + 11*x - 19;
f4 := 22*(x - 2)*y + (5*x^3 + 17*x^2 - 112*x + 120);
f5 := 11*y + (2*x^2 + 17*x - 34);
f6 := (x - 4)*y - (5*x - 9);

j := (f1*f2*f3*f4)^3/(f5^2*f6^11);

phi := map<X2 -> X1 | [j, 1]>;
return phi;