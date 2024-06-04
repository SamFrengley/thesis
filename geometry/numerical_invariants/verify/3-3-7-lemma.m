//////////////////////////////////////////////////
// Verify the statement about elements
// in GL_2(Z/2Z) and GL_2(Z/4Z) from Lemma 3.3.7
//////////////////////////////////////////////////

function InClasses(g, elts)
  G := Parent(g);
  all_classes := &join {Class(G, x) : x in elts};

  return g in all_classes;
end function;

for k in [1..2] do
  G := GL(2, Integers(2^k));

  I := [G![a,0,0,a] : a in [1..2^k] | a mod 2 eq 1];
  g_s := G![1,0,0,-1];
  g_ns := G![1,-2,2,-1];
  g_b := [G![a,2^(k-1),0,a] : a in [1..2^k] | a mod 2 eq 1];
  g_c := [G![0,-r,1,0] : r in [1..2^k] | r mod 2 eq 1];
  
  S := {g : g in G | g^2 in [G![Determinant(g),0,0,Determinant(g)], G![-Determinant(g),0,0,-Determinant(g)]]};

  for g in S do
    assert InClasses(g, I cat [g_s, g_ns] cat g_b cat g_c);
  end for;
end for;
