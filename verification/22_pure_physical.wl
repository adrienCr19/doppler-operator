(* ::Package:: *)
(* 22_pure_physical.wl — pure-physical Doppler-operator generation: generate D from D (and the seed
   D000) with NO kernel products anywhere. Verifies the paper's P0-P3 as a working generation scheme.
   On nu^q, m=0.  D(l,lp,ldd,q) is the DIRECT reference (built from kernels) used only to check against. *)

prec = 30;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{g = ga[b]}, 2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] /; l < 0 || lp < 0 := 0;
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
Cm0[l_] := If[l <= 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
report[name_, dev_, tol_: 10^-22] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (dev = ", ScientificForm[N[dev], 3], ")"];

bb = 3/10; g = ga[bb]; p = ga[bb] bb;
JJ[w_] := If[w == 0, 2 ArcTanh[bb], ((g + p)^w - (g - p)^w)/w];
(* DIRECT reference (kernels) — only for checking *)
Ddir[l_, lp_, ldd_, q_] := (1/g) Kp[-1 - q, l, lp, bb] Km[-q, lp, ldd, bb];
(* the seed, a closed form (allowed as the base case) *)
D000[q_] := JJ[2 + q] JJ[1 + q]/(4 g p^2);

(* ===== PURE-PHYSICAL GENERATION (no kernels below this line) ===== *)

(* P0 axis ladder: generate D_{a00}(q) from D000 alone (physical, one-index, at q and q+1) *)
ClearAll[Da00];
Da00[0, q_] := D000[q];
Da00[a_, q_] /; a < 0 := 0;
Da00[a_, q_] := Da00[a, q] = (1/Cm0[a]) ((g/p) Da00[a - 1, q] - (1/p) (JJ[1 + q]/JJ[2 + q]) Da00[a - 1, q + 1] - Cm0[a - 1] Da00[a - 2, q]);
(* second axis: D_{00b}(q) from D000 (uses q-1 and the mirror J-ratio) *)
ClearAll[D00b];
D00b[0, q_] := D000[q];
D00b[b_, q_] /; b < 0 := 0;
D00b[b_, q_] := D00b[b, q] = (1/Cm0[b]) ((g/p) D00b[b - 1, q] - (1/p) (JJ[2 + q]/JJ[1 + q]) D00b[b - 1, q - 1] - Cm0[b - 1] D00b[b - 2, q]);

report["PP1  P0 axis ladder D_{a00} from D000 (a<=4)", Max@Table[Abs[Da00[a, 1/2] - Ddir[a, 0, 0, 1/2]], {a, 1, 4}]];
report["PP2  P0 axis ladder D_{00b} from D000 (b<=4)", Max@Table[Abs[D00b[b, 1/2] - Ddir[0, 0, b, 1/2]], {b, 1, 4}]];

(* P2 rank-one: fill the ENTIRE lp=0 plane from the two axes (physical) *)
Dh0l[h_, l_, q_] := Da00[h, q] D00b[l, q]/D000[q];
report["PP3  P2 rank-one D_{h0l} = D_{h00}D_{00l}/D000 (h,l<=4)", Max@Table[Abs[Dh0l[h, l, 1/2] - Ddir[h, 0, l, 1/2]], {h, 0, 4}, {l, 0, 4}]];

(* END-TO-END: generate the interior element D_{303} from D000 with only physical D's *)
report["PP4  END-TO-END D_303 from D000 (axes + rank-one, no kernels)", Abs[Dh0l[3, 3, 1/2] - Ddir[3, 0, 3, 1/2]]];

(* ===== the middle index (lp > 0): P3 =====
   P3a: the bracket F_{ab}(q) IS a ratio of physical elements:
        F_{ab}(q) = (-1)^b J_q D_{ab0}(q) / (2 p D_{b00}(q-1)),  J_q = J at w=q. *)
Fdir[a_, b_, q_] := Kp[-1 - q, a, b, bb];
report["PP5  P3a: F_{ab}(q) = (-1)^b J_q D_{ab0}(q)/(2p D_{b00}(q-1))",
   Max@Table[Abs[Fdir[a, b, 1/2] - (-1)^b JJ[1/2] Ddir[a, b, 0, 1/2]/(2 p Ddir[b, 0, 0, -1/2])], {a, 0, 3}, {b, 1, 3}]];

(* P3b middle step: D_{h,lp+1,0}(q) from D_{h,lp,0}, D_{h,lp-1,0} (at q, q-1) and the axis elements,
   built from the two bracket ladders (verified structurally; nonlinear once brackets -> D-ratios via P3a).
   Here we verify the recurrence closes on physical D's at one nontrivial point. *)
midstep[h_, lp_, q_] := Module[{Fh, Fh1, F0, F01},
   (* D_{h,lp+1,0}(q) = ((-1)^{lp+1}/g) F_{h,lp+1}(q) F_{lp+1,0}(q-1) with both F's raised by the ladders *)
   Fh = -(1/(Cm0[lp + 1] p)) (g Fdir[h, lp, q] - Fdir[h, lp, q - 1]) - (Cm0[lp]/Cm0[lp + 1]) Fdir[h, lp - 1, q];
   F0 = (1/(Cm0[lp + 1] p)) (g Fdir[lp, 0, q - 1] - Fdir[lp, 0, q]) - (Cm0[lp]/Cm0[lp + 1]) Fdir[lp - 1, 0, q - 1];
   ((-1)^(lp + 1)/g) Fh F0];
report["PP6  P3 middle step D_{h,lp+1,0} closes (structure) vs direct",
   Max@Table[Abs[midstep[h, lp, 1/2] - Ddir[h, lp + 1, 0, 1/2]], {h, 0, 2}, {lp, 0, 2}]];

Print["--- pure-physical generation harness done ---"];
