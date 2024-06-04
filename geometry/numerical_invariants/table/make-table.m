AttachSpec("../spec");
SetAutoColumns(false);
SetColumns(1024);


table_str := "";
table_str := table_str cat "\\begin{tabular}{cc|ccccc} \n";
table_str := table_str cat "$N$ & $r$ & $p_g(\\ZNrtil{N}{r})$ & $\\kappa(\\ZNrtil{N}{r})$ & $p_g(\\WNr{N}{r})$ & $K_{W^{\\circ}}^2$ & $\\kappa(\\WNr{N}{r})$ \\\\ \n\\hline\n";

lb := 28;
ub := 40;

for N in [lb..ub] do
  rs := [r : r in [1..N] | GCD(r,N) eq 1];
  rr := [1];
  for r in rs do
    if not exists{rh : rh in rr | IsSquare(Integers(N)!(r*rh))} then
      Append(~rr, r);
    end if;
  end for;

  for r in rr do
    if r eq 1 then
      N_str := Sprintf("\\multirow{%o}{*}{$%o$}", #rr, N);
    else
      N_str := "";
    end if;
    pg_Z := GeometricGenusZNr(N,r); kd_Z := Min(2, pg_Z-1);
    pg_W := GeometricGenusWNr(N,r); kd_W := Min(2, pg_W-1);

    if pg_W ne 0 then
      K2 := Sprint(KWoSquared(N,r));
    else
      K2 := " ";
    end if;
    
    kd_W_str := Sprint(kd_W);
    
    if kd_W in [0,1] then
      if not N in [17,19] then
        kd_W_str := Sprintf("%o^{*}", kd_W);
      end if;

    elif pg_W ge 3 then
      if not KWoSquared(N,r) gt 0 then
        kd_W_str := Sprintf("%o^{*}", kd_W);
      end if;
    end if;

    table_str := table_str cat Sprintf("%o & $%o$ & $%o$ & $%o$ & $%o$ & $%o$ & $%o$ \\\\ \n", N_str ,r, pg_Z, kd_Z, pg_W, K2, kd_W_str);
    
  end for;

  if N ne ub then 
    table_str := table_str cat "\\hdashline \n";
  end if;
  
end for;

table_str := table_str cat "\\end{tabular}";

PrintFile("table_str.txt", table_str : Overwrite:=true);
