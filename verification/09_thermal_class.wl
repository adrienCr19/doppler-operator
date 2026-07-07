(* ::Package:: *)
(* 09_thermal_class.wl — Characterize the function class of thermally averaged Doppler operators.
   Claims:
   (L1) ladder: (d^2/dz^2 - 1) F(mu, a, z) = F(mu+1, a, z),  F(mu,a,z) = Int e^{-z cosh} cosh(a eta) sinh^{2mu}
        [immediate by cosh^2 = 1 + sinh^2; generalizes the sigma=2 telescoper of the companion work]
   (L2) sigma=2 channel (D101): G(z) := Sum_j c_j F(-1, a_j, z) satisfies
        (d^2/dz^2 - 1) G(z) = Sum_j c_j K_{a_j}(z)
   (L3) sigma=4 channel (D02):  (d^2/dz^2 - 1)^2 G(z) = Sum_j c_j K_{a_j}(z)
   So  K_2(z)/z * <D_{l l' l''}> is annihilated onto a finite Bessel combination by (d^2/dz^2-1)^{sigma/2},
   sigma = l + 2l' + l''.  Verified by high-precision finite differences in z. *)

prec = 40;
Ffun[mu_, a_, z_?NumericQ] := NIntegrate[Exp[-z Cosh[e]] Cosh[a e] Sinh[e]^(2 mu), {e, 0, Infinity},
   WorkingPrecision -> prec, PrecisionGoal -> 25, MaxRecursion -> 20];
report[name_, dev_, tol_] := Print[name, ": ", If[TrueQ[N[dev] < tol] || PossibleZeroQ[dev], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

(* channel data (Ohat -> -q), from 08_paper_tables.wl (verified against N8 closed forms) *)
tab101[q_] := With[{O = -q}, {
    {1, 3 (O^2 - 3 O + 6)/(8 O (O - 1) (O - 2) (O - 3))},
    {3, -3/(8 O (O - 3))},
    {1 - 2 O, -3/(8 O (O - 1))},
    {3 - 2 O, 3 (O^2 - 3 O + 1)/(4 O (O - 1) (O - 2) (O - 3))},
    {5 - 2 O, -3/(8 (O - 2) (O - 3))}}];
tab02[q_] := With[{O = -q, D1 = (-q - 1) (-q - 2)}, {
    {1, -5 (O^4 - 6 O^3 + 17 O^2 - 24 O + 72)/(16 O (O - 1)^2 (O - 2)^2 (O - 3))},
    {3, 15 (O^2 - 3 O + 6)/(32 D1^2)},
    {5, -5/(32 D1)},
    {1 - 2 O, -5 (2 O^3 - 12 O^2 + 13 O + 12)/(16 O D1^2)},
    {3 - 2 O, 15 (O^2 - 3 O - 3)/(16 D1^2)},
    {5 - 2 O, -5 (O + 1) (2 O^2 - 8 O + 3)/(16 (O - 3) D1^2)},
    {7 - 2 O, 5 O (O + 1)/(32 (O - 1) (O - 2)^2 (O - 3))},
    {1 + 2 O, 5 (O - 3) (O - 4)/(32 O (O - 1)^2 (O - 2))}}];

(* combined (regularized) channel integrals *)
Gchan[tab_, sig_, z_?NumericQ] := NIntegrate[
   Exp[-z Cosh[e]] (Total[tab /. {a_?NumericQ, c_} :> c Cosh[a e]])/Sinh[e]^sig,
   {e, 0, Infinity}, WorkingPrecision -> prec, PrecisionGoal -> 25, MaxRecursion -> 20];

(* high-order central second derivative, 7-point, h chosen for 40-digit arithmetic *)
d2[f_, z_, h_] := (2 f[z - 3 h] - 27 f[z - 2 h] + 270 f[z - h] - 490 f[z] + 270 f[z + h] - 27 f[z + 2 h] + 2 f[z + 3 h])/(180 h^2);

z0 = 20; h = 1/100;  (* z = 1/theta = 20  <->  kT ~ 25.6 keV *)
qv = 1/2;

(* ---- L1: base of the ladder. NOTE: individual F(mu<0, a, z) diverge at eta -> 0 (csch pole);
   only the Sum_j c_j = 0 regularized channel combinations exist, tested in L2/L3. ---- *)
dev = Max@Table[Abs[Ffun[0, a, z0] - BesselK[a, z0]], {a, {1, 3 + 2 qv}}];
report["L1  F(0,a,z) == BesselK (sigma = 0 rung)", dev, 10^-25];

(* ---- L2: sigma=2 channel annihilated to Bessel in one step ---- *)
bess[tab_, z_] := Total[tab /. {a_?NumericQ, c_} :> c BesselK[a, z]];
G101[z_?NumericQ] := Gchan[tab101[qv], 2, z];
dev = Abs[d2[G101, z0, h] - G101[z0] - bess[tab101[qv], z0]];
report["L2  (d2-1) G_101(z) = Sum_j c_j K_{a_j}(z)", dev, 10^-15];

(* ---- L3: sigma=4 channel annihilated in two steps ---- *)
G02[z_?NumericQ] := Gchan[tab02[qv], 4, z];
step1[z_?NumericQ] := d2[G02, z, h] - G02[z];              (* = Sum c_j F(-1, a_j, z) *)
dev1 = Abs[step1[z0] - Total[tab02[qv] /. {a_?NumericQ, c_} :> c Ffun[-1, a, z0]]];
report["L3a (d2-1) G_02 = Sum_j c_j F(-1,a_j,z)", dev1, 10^-13];
dev2 = Abs[d2[step1, z0, 1/50] - step1[z0] - bess[tab02[qv], z0]];
report["L3b (d2-1)^2 G_02(z) = Sum_j c_j K_{a_j}(z)", dev2, 10^-9];

Print["--- thermal-class harness done ---"];
