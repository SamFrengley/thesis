/*

*/
Attach("../ZNr-equations.m");
Attach("/home/stf32/Documents/MagmaFiles/Z99_OthersCode/RSZB_ladic/groups/gl2.m");
AttachSpec("/mhome/dpmms/t/stf32/Documents/MagmaFiles/Z99_OthersCode/Zywina_OpenImage/OpenImage.spec");
Z := OpenImageContext("/mhome/dpmms/t/stf32/Documents/MagmaFiles/Z99_OthersCode/Zywina_OpenImage/data-files");

//---------------------------
// Isogenous Claim
t := 2;
xi, eta := Explode([0,4]);

aa := eval Read("../Z12-r/infinitefamilies/12-1.m");
E1_12_1 := EllipticCurve(aa[1]);
assert not jInvariant(E1_12_1) in {0, 1728};
E1_12_11 := MinimalTwist(E1_12_1);
assert #IsogenousCurves(E1_12_1) eq 1;

aa := eval Read("../Z12-r/infinitefamilies/12-5.m");
E1_12_5 := EllipticCurve(aa[1]);
assert not jInvariant(E1_12_5) in {0, 1728};
E1_12_5 := MinimalTwist(E1_12_5);
assert #IsogenousCurves(E1_12_5) eq 1;

aa := eval Read("../Z12-r/infinitefamilies/12-7.m");
E1_12_7 := EllipticCurve(aa[1]);
assert not jInvariant(E1_12_7) in {0, 1728};
E1_12_11 := MinimalTwist(E1_12_7);
assert #IsogenousCurves(E1_12_7) eq 1;

aa := eval Read("../Z12-r/infinitefamilies/12-11.m");
E1_12_11 := EllipticCurve(aa[1]);
assert not jInvariant(E1_12_11) in {0, 1728};
E1_12_11 := MinimalTwist(E1_12_11);
assert #IsogenousCurves(E1_12_11) eq 1;

aa := eval Read("../Z14-r/infinitefamilies/14-1.m");
E1_14_1 := EllipticCurve(aa[1]);
assert not jInvariant(E1_14_1) in {0, 1728};
E1_14_1 := MinimalTwist(E1_14_1);
assert #IsogenousCurves(E1_14_1) eq 1;


//---------------------------
// Surjectivity Claim
function IsSurjectiveGaloisRep(E, N)
  img := FindOpenImage(Z, E);
  M := Characteristic(BaseRing(img));
  modN_img := GL2Project(img, GCD(M,N));
  if #modN_img eq #GL(2, Integers(GCD(M,N))) then
    return true;
  else
    return false;
  end if;
end function;

assert IsSurjectiveGaloisRep(E1_12_1, 12);
assert IsSurjectiveGaloisRep(E1_12_5, 12);
assert IsSurjectiveGaloisRep(E1_12_7, 12);
assert IsSurjectiveGaloisRep(E1_12_11, 12);

assert IsSurjectiveGaloisRep(E1_14_1, 14);
