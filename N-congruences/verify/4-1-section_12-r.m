/*
    Verifying the Hecke correspondences in Tables 5.2--5.9
*/
load "../Z12-r/Z12equations.m";
QQ := Rationals();
A3<u,v,z> := AffineSpace(QQ, 3);

function CorrectMap(m, ff)
// This returns the map from X_0(m) (in the small modular curve database) to the copy of 
// X_0^+(m) on our model for Z(12, m). This will be used to check the curves we have 
// do map to Hecke correspondences
    if m eq 5 then
        x0 := Explode(ff);
        pp := [
                0,
                16*x0/(x0^2 + 22*x0 + 125)
        ];
        return pp;
    
    elif m eq 7 then
        x0 := Explode(ff);
        pp := [
                0,
                (-1/2*x0^2 - 13/2*x0 - 49/2)/x0
            ];
        return pp;
    
    elif m eq 11 then
        x0, x1 := Explode(ff);
        pp := [
                4/(x0^2 + 6*x0 + 11)*x1 + (x0^2 - 6*x0 - 37)/(x0^2 + 6*x0 + 11),
                -1
            ];
        return pp;

    elif m eq 13 then
        x0 := ff[1];
        pp := [
                -2*x0/(x0^2 + 6*x0 + 13),
                -1
            ];
        return pp;
    
    elif m eq 17 then
        x0, x1 := Explode(ff);
        pp := [
                1/(x0^2 + 4*x0 + 8)*x1 + (x0^2 + 3*x0 + 2)/(x0^2 + 4*x0 + 8),
                -1/(x0^2 + 4*x0 + 8)*x1 + (-x0^2 - 3*x0 - 2)/(x0^2 + 4*x0 + 8)
            ];
        return pp;
    
    elif m eq 19 then
        x0, x1 := Explode(ff);
        pp := [
                -2/(x0^2 + 2*x0 + 3)*x1 + (4*x0 - 2)/(x0^2 + 2*x0 + 3),
                -1
            ];
        return pp;

    elif m eq 23 then
        x0, x1 := Explode(ff);
        pp := [
                (x0^2 + 3)/(x0^2 - 1),
                (-x0 + 1)/(x0 + 1)
            ];
        return pp;

    elif m eq 25 then
        x0 := ff[1];
        pp := [
                -2*x0^2/(x0^4 + 4*x0^3 + 14*x0^2 + 20*x0 + 25),
                (-x0^2 - 4*x0 - 5)/(x0^2 + 2*x0 + 5)
            ];
        return pp;

    elif m eq 29 then
        x0, x1 := Explode(ff);
        pp := [
                (x0^2 - x0 + 1)/(x0^2 - x0),
                (x0^2 + x0 - 1)/(x0^2 - x0)
            ];
        return pp;

    elif m eq 31 then
        x0, x1 := Explode(ff);
        pp := [
                -2/(x0 - 1),
                (x0^2 - x0 - 2)/(x0^2 - x0)
            ];
        return pp;

    elif m eq 35 then
        x0, x1 := Explode(ff);
        pp := [
                (x0^4 + x0^3 + 3*x0 - 1)/(x0^4 + x0^3 - x0 - 1),
                (x0 - 1)/(x0 + 1)
            ];
        return pp;

    elif m eq 37 then
        x0, x1 := Explode(ff);
        pp := [
            (-2*x0^3 - 2*x0^2 + 2*x0)/(20*x0^6 - 36*x0^5 + 57*x0^4 - 46*x0^3 + 33*x0^2 -
                12*x0 + 4)*x1 + (2*x0^6 + 6*x0^5 - 16*x0^4 + 24*x0^3 - 22*x0^2 + 12*x0 -
                4)/(20*x0^6 - 36*x0^5 + 57*x0^4 - 46*x0^3 + 33*x0^2 - 12*x0 + 4),
            6*x0/(10*x0^4 - 13*x0^3 + 17*x0^2 - 8*x0 + 4)*x1 + (4*x0^4 - 5*x0^3 + 
                5*x0^2)/(10*x0^4 - 13*x0^3 + 17*x0^2 - 8*x0 + 4)
            ];
        return pp;

    elif m eq 41 then
        x0, x1 := Explode(ff);
        pp := [
                (x0^3 - 2*x0^2 + 2)/(x0^3 - 2*x0^2 + 1),
                (x0^3 - 2*x0)/(x0^3 - 2*x0^2 + 1)
            ];
        return pp;

    elif m eq 43 then
        x0, x1 := Explode(ff);
        pp := [
                (2*x0^2 - 8*x0 + 6)/(x0^4 - x0^3 + 7*x0^2 - 2*x0 + 10)*x1^2 + (2*x0^3 
                - 14*x0^2 + 12*x0 - 8)/(x0^4 - x0^3 + 7*x0^2 - 2*x0 + 10)*x1 + (-10*x0^3 
                + 6*x0^2 - 18*x0 + 6)/(x0^4 - x0^3 + 7*x0^2 - 2*x0 + 10),
                (2*x0 - 2)/(x0^2 - x0 + 5)*x1^2 + (4*x0^2 - 6*x0 + 6)/(x0^2 - x0 + 5)*x1 
                + (2*x0^3 - 5*x0^2 + 3*x0 - 7)/(x0^2 - x0 + 5)
            ];
        return pp;

    elif m eq 47 then
        x0, x1 := Explode(ff);
        pp := [
            (x0^4 + x0^2 + 2)/(x0^4 + x0^2 - 2),
            (-x0 + 1)/(x0 + 1)
        ];
        return pp;

    elif m eq 49 then
        x0, x1 := Explode(ff);
        pp := [
            (8*x0^4 - 24*x0^3 + 8*x0^2 + 32)/(x0^8 - 4*x0^7 + 16*x0^6 - 4*x0^5 + 6*x0^4 
                + 100*x0^3 + 160*x0^2 + 100*x0 + 25)*x1 + (-2*x0^6 + 28*x0^4 - 32*x0^3 -
                42*x0^2 + 40*x0 + 24)/(x0^8 - 4*x0^7 + 16*x0^6 - 4*x0^5 + 6*x0^4 + 
                100*x0^3 + 160*x0^2 + 100*x0 + 25),
            (-2*x0^2 + 4*x0 - 14)/(x0^4 - 2*x0^3 + 6*x0^2 + 10*x0 + 5)*x1 + (-x0^4 + 
                2*x0^3 - 4*x0^2 - 22*x0 - 3)/(x0^4 - 2*x0^3 + 6*x0^2 + 10*x0 + 5)
        ];
        return pp;

    elif m eq 53 then
        x0, x1, x2 := Explode(ff);
        pp := [ 
            (x2 - 1)/(2*x2^4 - 3*x2^3 + 5*x2^2 - 3*x2 + 2)*x1^3 + (x2^2 - 2*x2)/(2*x2^4 
                - 3*x2^3 + 5*x2^2 - 3*x2 + 2)*x1^2 + (x2^3 - 2*x2^2 + x2 - 1)/(2*x2^4 - 
                3*x2^3 + 5*x2^2 - 3*x2 + 2)*x1 + (2*x2^5 - x2^4 - 2*x2^3 + 7*x2^2 - 5*x2
                + 2)/(2*x2^5 - 3*x2^4 + 5*x2^3 - 3*x2^2 + 2*x2),
            (x2 - 1)/(2*x2^4 - 3*x2^3 + 5*x2^2 - 3*x2 + 2)*x1^3 + (x2^2 - 2*x2)/(2*x2^4 
                - 3*x2^3 + 5*x2^2 - 3*x2 + 2)*x1^2 + (x2^3 - 2*x2^2 + x2 - 1)/(2*x2^4 - 
                3*x2^3 + 5*x2^2 - 3*x2 + 2)*x1 + (2*x2^5 - 5*x2^4 + 4*x2^3 - 3*x2^2 + x2
                - 2)/(2*x2^5 - 3*x2^4 + 5*x2^3 - 3*x2^2 + 2*x2)
        ];
        return pp;

    elif m eq 55 then
        x0, x1 := Explode(ff);
        pp := [
            (-2*x1^3 - 2*x1^2)/(15*x1^4 - 10*x1^3 - 35*x1^2 + 20*x1 + 10)*x0^3 + (4*x1^4
                + 2*x1^3 - 4*x1^2 - 4*x1 - 2)/(15*x1^4 - 10*x1^3 - 35*x1^2 + 20*x1 + 
                10)*x0^2 + (-2*x1^5 - 16*x1^4 + 2*x1^3 + 20*x1^2 - 2*x1 - 6)/(15*x1^4 - 
                10*x1^3 - 35*x1^2 + 20*x1 + 10)*x0 + (-2*x1^5 - 10*x1^4 + 18*x1^3 + 
                12*x1^2 + 2*x1 - 4)/(15*x1^4 - 10*x1^3 - 35*x1^2 + 20*x1 + 10),
            (x1 + 1)/(x1 - 1)
        ];
        return pp;

    elif m eq 59 then
        x0, x1 := Explode(ff);
        pp := [
                (x0^5 - x0^4 + x0^2 - x0 + 4)/(x0^5 - x0^4 + x0^2 - x0),
                (-x0 + 1)/(x0 + 1)
        ];
        return pp;

    elif m eq 61 then
        x0, x1, x2 := Explode(ff);
        pp := [
            (2*x1^4 - 4*x1^2 - 2*x1)/(2*x1^7 + 13*x1^6 + 30*x1^5 + 36*x1^4 + 36*x1^3 + 
                24*x1^2 + 8*x1 + 1)*x2^3 + (-2*x1^5 + 6*x1^3 + 2*x1^2 - 4*x1 - 
                2)/(2*x1^7 + 13*x1^6 + 30*x1^5 + 36*x1^4 + 36*x1^3 + 24*x1^2 + 8*x1 + 
                1)*x2^2 + (2*x1^6 + 6*x1^5 - 2*x1^4 - 16*x1^3 - 10*x1^2 + 2*x1 + 
                2)/(2*x1^7 + 13*x1^6 + 30*x1^5 + 36*x1^4 + 36*x1^3 + 24*x1^2 + 8*x1 + 
                1)*x2 + (-2*x1^6 - 2*x1^5 - 6*x1^4 - 2*x1^3 + 2*x1^2 + 2*x1)/(x1^6 + 
                6*x1^5 + 12*x1^4 + 12*x1^3 + 12*x1^2 + 6*x1 + 1),
            (2*x1^2 + 2*x1)/(2*x1^5 + 7*x1^4 + 7*x1^3 + 8*x1^2 + 5*x1 + 1)*x2^3 + 
                (-2*x1^3 - 2*x1^2 + 2*x1 + 2)/(2*x1^5 + 7*x1^4 + 7*x1^3 + 8*x1^2 + 5*x1 
                + 1)*x2^2 + (2*x1^4 + 8*x1^3 + 8*x1^2 - 2)/(2*x1^5 + 7*x1^4 + 7*x1^3 + 
                8*x1^2 + 5*x1 + 1)*x2 + (-x1^4 + x1^3 - 2*x1^2 - x1 + 1)/(x1^4 + 3*x1^3 
                + 2*x1^2 + 3*x1 + 1)
        ];
        return pp;

    elif m eq 71 then
        x0, x1 := Explode(ff);
        pp := [
            (x0^6 + 4*x0^5 + 5*x0^4 - 5*x0^2 - 2*x0 + 4)/(x0^6 + 4*x0^5 + 5*x0^4 - 
                5*x0^2 - 2*x0),
            -x0/(x0 + 2)
        ];
        return pp;

    end if;

    return 0;
end function;


// -------------------------------------
// We now check that the modular curves we claim actually map to 
// Hecke correspondences on Z(1) do so
// -------------------------------------

//---------------------
// (12, 1)
ZZ := Surface(A3, z^2 - Z12Equations(1, u, v));
table2 := eval Read("../Z12-r/curvesonZ12-r/modularcurves/Z12-1.m");

for ex in table2 do
    m, X := Explode(ex);
    printf "The case when m=%o\n", m;
    X := Curve(ZZ, [X]);
    X0m := SmallModularCurve(m);
    g := Genus(X0m);
    assert Genus(X) eq g;

    if g eq 0 then
        KK := FunctionField(QQ);
    else
        KK := FunctionField(X0m);
    end if;

    X0m := BaseChange(X0m, KK);
    coords := [KK.i : i in [1..Dimension(AmbientSpace(X0m))]];
    gen_pt := coords cat [1];

    pp := CorrectMap(m, coords);

    jj, jj_1728 := W12Moduli(1, pp[1],pp[2]);
    j1 := jNInvariant(X0m!gen_pt, m);
    j2 := jInvariant(X0m!gen_pt, m);

    assert jj eq j1*j2;
    assert jj_1728 eq (j1 - 1728)*(j2 - 1728);
end for;


//---------------------
// (12, 5)
ZZ := Surface(A3, z^2 - Z12Equations(5, u, v));

// The blown down case when m=5
printf "The case when m=5\n";
Kt<t> := FunctionField(QQ);
P<e> := PowerSeriesRing(Kt);
A2 := AffineSpace(QQ, 2);

bl_up := Z12Equations(5, 1/2*(3 - e), 1/2*(1 + e + e^2*t));
assert bl_up + BigO(e^5) eq -1/4*(t^2 + 44*t - 16)*e^4 + BigO(e^5);

jj, jj_1728 := W12Moduli(5, 1/2*(3 - e), 1/2*(1 + e + e^2*t));
jj := Evaluate(jj,0);
jj_1728 := Evaluate(jj_1728,0);

X0m := SmallModularCurve(5);
X0m := BaseChange(X0m, Kt);
pp := CorrectMap(5, [t]);
assert pp[1] eq 0;

jj := Evaluate(jj, pp[2]);
jj_1728 := Evaluate(jj_1728, pp[2]);
j1 := jNInvariant(X0m![t,1], 5);
j2 := jInvariant(X0m![t,1], 5);

assert jj eq j1*j2;
assert jj_1728 eq (j1 - 1728)*(j2 - 1728);

// The remaining cases when r=5
table3 := eval Read("../Z12-r/curvesonZ12-r/modularcurves/Z12-5.m");

for ex in table3 do
    m, X := Explode(ex);
    printf "The case when m=%o\n", m;
    X := Curve(ZZ, [X]);
    X0m := SmallModularCurve(m);
    g := Genus(X0m);
    assert Genus(X) eq g;

    if g eq 0 then
        KK := FunctionField(QQ);
    else
        KK := FunctionField(X0m);
    end if;

    X0m := BaseChange(X0m, KK);
    coords := [KK.i : i in [1..Dimension(AmbientSpace(X0m))]];
    gen_pt := coords cat [1];

    pp := CorrectMap(m, coords);

    jj, jj_1728 := W12Moduli(5, pp[1],pp[2]);
    j1 := jNInvariant(X0m!gen_pt, m);
    j2 := jInvariant(X0m!gen_pt, m);

    assert jj eq j1*j2;
    assert jj_1728 eq (j1 - 1728)*(j2 - 1728);
end for;


//---------------------
// (12, 7)
ZZ := Surface(A3, z^2 - Z12Equations(7, u, v));

// The blown down case when m=7
printf "The case when m=7\n";
Kt<t> := FunctionField(QQ);
P<e> := PowerSeriesRing(Kt);
A2 := AffineSpace(QQ, 2);

bl_up := Z12Equations(7, -e, 1 + e + e^2*t);
assert bl_up + BigO(e^5) eq (4*t^2 + 52*t - 27)*e^4 + BigO(e^5);

jj, jj_1728 := W12Moduli(7, -e, 1 + e + e^2*t);
jj := Evaluate(jj,0);
jj_1728 := Evaluate(jj_1728,0);

X0m := SmallModularCurve(7);
X0m := BaseChange(X0m, Kt);
pp := CorrectMap(7, [t]);
assert pp[1] eq 0;

jj := Evaluate(jj, pp[2]);
jj_1728 := Evaluate(jj_1728, pp[2]);
j1 := jNInvariant(X0m![t,1], 7);
j2 := jInvariant(X0m![t,1], 7);

assert jj eq j1*j2;
assert jj_1728 eq (j1 - 1728)*(j2 - 1728);

// The remaining cases when r=7
table4 := eval Read("../Z12-r/curvesonZ12-r/modularcurves/Z12-7.m");

for ex in table4 do
    m, X := Explode(ex);
    printf "The case when m=%o\n", m;
    X := Curve(ZZ, [X]);
    X0m := SmallModularCurve(m);
    g := Genus(X0m);
    assert Genus(X) eq g;

    if g eq 0 then
        KK := FunctionField(QQ);
    else
        KK := FunctionField(X0m);
    end if;

    X0m := BaseChange(X0m, KK);
    coords := [KK.i : i in [1..Dimension(AmbientSpace(X0m))]];
    gen_pt := coords cat [1];

    pp := CorrectMap(m, coords);

    jj, jj_1728 := W12Moduli(7, pp[1],pp[2]);
    j1 := jNInvariant(X0m!gen_pt, m);
    j2 := jInvariant(X0m!gen_pt, m);

    assert jj eq j1*j2;
    assert jj_1728 eq (j1 - 1728)*(j2 - 1728);
end for;


//---------------------
// (12, 11)
ZZ := Surface(A3, z^2 - Z12Equations(11, u, v));

table5 := eval Read("../Z12-r/curvesonZ12-r/modularcurves/Z12-11.m");

for ex in table5 do
    m, X := Explode(ex);
    printf "The case when m=%o\n", m;
    X := Curve(ZZ, [X]);
    X0m := SmallModularCurve(m);
    g := Genus(X0m);
    assert Genus(X) eq g;

    if g eq 0 then
        KK := FunctionField(QQ);
    else
        KK := FunctionField(X0m);
    end if;

    X0m := BaseChange(X0m, KK);
    coords := [KK.i : i in [1..Dimension(AmbientSpace(X0m))]];
    gen_pt := coords cat [1];

    pp := CorrectMap(m, coords);

    jj, jj_1728 := W12Moduli(11, pp[1],pp[2]);
    j1 := jNInvariant(X0m!gen_pt, m);
    j2 := jInvariant(X0m!gen_pt, m);

    assert jj eq j1*j2;
    assert jj_1728 eq (j1 - 1728)*(j2 - 1728);
end for;
