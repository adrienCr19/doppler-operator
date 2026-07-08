(* ::Package:: *)
(* 11_f_calculus.wl — The F-calculus (N11): complete recurrence toolkit for the thermal-average
   special functions
     F(mu,a,z) = Int_0^inf e^{-z cosh eta} cosh(a eta) sinh^{2 mu} eta d eta   (even sector)
     L(mu,a,z) = Int_0^inf e^{-z cosh eta} sinh(a eta) sinh^{2 mu} eta d eta   (odd sector)
   Relations (derived on proofs/f-calculus.html):
     FC1  (d2/dz2 - 1) F(mu,a,z) = F(mu+1,a,z);  same for L                       [mu-ladder]
     FC2  F(mu,a+1) + F(mu,a-1) = -2 dF(mu,a)/dz;  same for L                     [a-sum]
     FC3  F(mu,a+1) - F(mu,a-1) = 2 L(mu+1/2, a);  L(mu,a+1) - L(mu,a-1) = 2 F(mu+1/2, a)   [a-difference]
     FC4  z F(mu+1,a) = a L(mu+1/2,a) - (2mu+1) dF(mu,a)/dz     (mu > -1/2)
     FC4b z F(1/2,a)  = a L(0,a) + e^{-z}                        (boundary-corrected mu = -1/2 case)
     FC4' z L(mu+1,a) = a F(mu+1/2,a) - (2mu+1) dL(mu,a)/dz     (mu >= -1/2; no boundary term)
     FC5  special values: F(0,a) = K_a(z);  L(1/2,a) = (a/z) K_a(z);
          F(1,a) = (a^2/z^2) K_a - K_a'/z    [ladder + Bessel ODE]
     FC6  channel-regularized FC4 at mu = -1 on the D101 channel data. *)

prec = 40;
FF[mu_, a_, z_?NumericQ] := NIntegrate[Exp[-z Cosh[e]] Cosh[a e] Sinh[e]^(2 mu), {e, 0, Infinity},
   WorkingPrecision -> prec, PrecisionGoal -> 22, MaxRecursion -> 20];
LL[mu_, a_, z_?NumericQ] := NIntegrate[Exp[-z Cosh[e]] Sinh[a e] Sinh[e]^(2 mu), {e, 0, Infinity},
   WorkingPrecision -> prec, PrecisionGoal -> 22, MaxRecursion -> 20];
d1[f_, z_, h_] := (-f[z + 2 h] + 8 f[z + h] - 8 f[z - h] + f[z - 2 h])/(12 h);
d2[f_, z_, h_] := (2 f[z - 3 h] - 27 f[z - 2 h] + 270 f[z - h] - 490 f[z] + 270 f[z + h] - 27 f[z + 2 h] + 2 f[z + 3 h])/(180 h^2);
report[name_, dev_, tol_: 10^-14] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

z0 = 12; h = 1/100; as = {3/2, 3 - 2 (1/2)};  (* generic non-integer orders incl. a = integer - 2*Ohat type *)

(* FC1 mu-ladder for F and L *)
dev = Max@Table[Abs[d2[FF[mu, a, #] &, z0, h] - FF[mu, a, z0] - FF[mu + 1, a, z0]], {mu, {0, 1/2, 1}}, {a, as}];
report["FC1a (d2-1)F(mu) = F(mu+1)", dev];
dev = Max@Table[Abs[d2[LL[mu, a, #] &, z0, h] - LL[mu, a, z0] - LL[mu + 1, a, z0]], {mu, {0, 1/2}}, {a, as}];
report["FC1b (d2-1)L(mu) = L(mu+1)", dev];

(* FC2 a-sum *)
dev = Max@Table[Abs[FF[mu, a + 1, z0] + FF[mu, a - 1, z0] + 2 d1[FF[mu, a, #] &, z0, h]], {mu, {0, 1/2, 1}}, {a, as}];
report["FC2a F(mu,a+1)+F(mu,a-1) = -2 F'(mu,a)", dev];
dev = Max@Table[Abs[LL[mu, a + 1, z0] + LL[mu, a - 1, z0] + 2 d1[LL[mu, a, #] &, z0, h]], {mu, {0, 1/2}}, {a, as}];
report["FC2b L version", dev];

(* FC3 a-difference (couples F and L, half-step in mu) *)
dev = Max@Table[Abs[FF[mu, a + 1, z0] - FF[mu, a - 1, z0] - 2 LL[mu + 1/2, a, z0]], {mu, {0, 1/2, 1}}, {a, as}];
report["FC3a F(mu,a+1)-F(mu,a-1) = 2 L(mu+1/2,a)", dev, 10^-20];
dev = Max@Table[Abs[LL[mu, a + 1, z0] - LL[mu, a - 1, z0] - 2 FF[mu + 1/2, a, z0]], {mu, {0, 1/2}}, {a, as}];
report["FC3b L(mu,a+1)-L(mu,a-1) = 2 F(mu+1/2,a)", dev, 10^-20];

(* FC4 z-relation, generic mu *)
dev = Max@Table[Abs[z0 FF[mu + 1, a, z0] - a LL[mu + 1/2, a, z0] + (2 mu + 1) d1[FF[mu, a, #] &, z0, h]], {mu, {0, 1/2, 1}}, {a, as}];
report["FC4  z F(mu+1) = a L(mu+1/2) - (2mu+1) F'(mu)", dev];
(* FC4b boundary-corrected mu = -1/2 *)
dev = Max@Table[Abs[z0 FF[1/2, a, z0] - a LL[0, a, z0] - Exp[-z0]], {a, as}];
report["FC4b z F(1/2,a) = a L(0,a) + e^{-z}", dev, 10^-20];
(* FC4' odd version incl. mu = -1/2 (no boundary term) *)
dev = Max@Table[Abs[z0 LL[mu + 1, a, z0] - a FF[mu + 1/2, a, z0] + (2 mu + 1) d1[LL[mu, a, #] &, z0, h]], {mu, {0, 1/2}}, {a, as}];
report["FC4' z L(mu+1) = a F(mu+1/2) - (2mu+1) L'(mu)", dev];
dev = Max@Table[Abs[z0 LL[1/2, a, z0] - a FF[0, a, z0]], {a, as}];
report["FC4'b z L(1/2,a) = a K_a(z)   [special value]", dev, 10^-20];

(* FC5 special values and first closed rungs *)
dev = Max@Table[Abs[FF[0, a, z0] - BesselK[a, z0]], {a, as}];
report["FC5a F(0,a) = K_a", dev, 10^-20];
dev = Max@Table[Abs[FF[1, a, z0] - ((a^2/z0^2) BesselK[a, z0] + (BesselK[a - 1, z0] + BesselK[a + 1, z0])/(2 z0))], {a, as}];
report["FC5b F(1,a) = (a^2/z^2)K_a - K_a'/z  (K_a' = -(K_{a-1}+K_{a+1})/2)", dev, 10^-20];
(* F(3/2) from the ladder of the elementary F(1/2): F(3/2) = (d2-1)F(1/2), F(1/2) = (a L(0,a)+e^-z)/z *)
F12[a_, z_?NumericQ] := (a LL[0, a, z] + Exp[-z])/z;
dev = Max@Table[Abs[FF[3/2, a, z0] - (d2[F12[a, #] &, z0, h] - F12[a, z0])], {a, as}];
report["FC5c F(3/2,a) = (d2-1)[(a L_a + e^-z)/z]  (ladder from elementary rung)", dev];

(* FC6 channel-regularized FC4 at mu = -1 (sigma = 2, D101 channel):
   z Sum_j c_j F(0,a_j) = Sum_j c_j a_j L(-1/2,a_j) + d/dz Sum_j c_j F(-1,a_j)
   [L(-1/2,a) is individually convergent; Sum c_j F(-1,a_j) evaluated as combined integrand] *)
qv = 1/2; Ov = -qv;
tab101 = {{1, 3 (Ov^2 - 3 Ov + 6)/(8 Ov (Ov - 1) (Ov - 2) (Ov - 3))}, {3, -3/(8 Ov (Ov - 3))},
   {1 - 2 Ov, -3/(8 Ov (Ov - 1))}, {3 - 2 Ov, 3 (Ov^2 - 3 Ov + 1)/(4 Ov (Ov - 1) (Ov - 2) (Ov - 3))},
   {5 - 2 Ov, -3/(8 (Ov - 2) (Ov - 3))}};
Gm1[z_?NumericQ] := NIntegrate[Exp[-z Cosh[e]] (Total[tab101 /. {a_?NumericQ, c_} :> c Cosh[a e]])/Sinh[e]^2,
   {e, 0, Infinity}, WorkingPrecision -> prec, PrecisionGoal -> 22];
lhs = z0 Total[tab101 /. {a_?NumericQ, c_} :> c BesselK[a, z0]];
rhs = Total[tab101 /. {a_?NumericQ, c_} :> c a LL[-1/2, a, z0]] + d1[Gm1, z0, h];
report["FC6  channel-regularized FC4 at mu = -1 (D101 data)", Abs[lhs - rhs]];

Print["--- F-calculus harness done ---"];
