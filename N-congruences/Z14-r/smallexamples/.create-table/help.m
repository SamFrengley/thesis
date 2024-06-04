SetAutoColumns(false);
SetColumns(1024);

congs := eval Read("14-all.m");
congs := [x : x in congs | x[2] lt 10^10 and x[3] lt 10^10];
table := [];
jj := {};
for x in congs do
  j := {jInvariant(EllipticCurve(x[5][1])),jInvariant(EllipticCurve(x[5][2]))};
  if not j in jj then
    Include(~jj, j);
    Append(~table, x);
  end if;
end for;
           
PrintFile("table.m", congs : Overwrite:=true);

if #table mod 2 eq 1 then
    table := &cat [[table[i], table[ExactQuotient(#table+1, 2) + i]] : i in [1..ExactQuotient(#table-1, 2)]] cat [table[ExactQuotient(#table+1, 2)]];
else
    table := &cat [[table[i], table[ExactQuotient(#table, 2) + i]] : i in [1..ExactQuotient(#table, 2)]];
end if;

str := "";

for i in [1..#table] do
    ex := table[i];
    r := ex[1];
    E1, E2 := Explode(ex[5]);
    E1 := EllipticCurve(E1); E2 := EllipticCurve(E2);

    if Conductor(E1) le 500000 then
        s1 := "\\LMFDBLabel{" cat LMFDBLabel(E1) cat "}";
    else
        s1 := "\\texttt{" cat Sprint(Conductor(E1)) cat "*}";
    end if;

    if Conductor(E2) le 500000 then
        s2 := "\\LMFDBLabel{" cat LMFDBLabel(E2) cat "}";
    else
        s2 := "\\texttt{" cat Sprint(Conductor(E2)) cat "*}";
    end if;

    new_st := Sprintf("%-4o & %-30o & %-30o ", "$" cat Sprint(r) cat "$", s1, s2);
    str := str cat new_st;

    if i mod 2 eq 1 then
        str := str cat " & ";
    else
        str := str cat "\\\\ \n";
    end if;
end for;

if #table mod 2 eq 1 then
    str := str cat Sprintf("%-4o & %-30o & %-30o ", "", "", "") cat "\\\\";
end if;

PrintFile(".me_tex.txt", str, "Minimal" : Overwrite:=true);
