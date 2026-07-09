(* ::Package:: *)
(* 23_crosscheck.wl — do the pure-physical generation results (harness 22 / paper P0-P4) agree with the
   previously-published CLOSED FORMS (N8 harness 06, N16 harness 16)?  Generate the flagship elements by the
   pure-physical scheme and compare to the closed forms. m=0, on nu^q. *)

prec = 30;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{gg = ga[b]}, 2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(gg (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] /; l < 0 || lp < 0 := 0;
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
Cm0[l_] := If[l <= 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
(* evaluate at 30 digits: C2/C3 compare two EXACT algebraic expressions (P0 ladder is built from the exact
   seed, no quadrature) — machine-precision N would show spurious ~1e-16 cancellation, so force high precision;
   quadrature-based devs (C4-C6) are unaffected. *)
report[name_, dev_, tol_: 10^-22] := Print[name, ": ", If[TrueQ[N[dev, 30] < tol], "PASS", "FAIL"], "  (dev = ", ScientificForm[N[dev, 30], 3], ")"];

bb = 3/10; g = ga[bb]; p = ga[bb] bb;
JJ[w_] := If[w == 0, 2 ArcTanh[bb], ((g + p)^w - (g - p)^w)/w];
Fd[a_, b_, q_] := Kp[-1 - q, a, b, bb];

(* ---------- PUBLISHED CLOSED FORMS (from N8 harness 06, N16 harness 16) ---------- *)
D000c[q_] := JJ[2 + q] JJ[1 + q]/(4 g p^2);
D101c[q_] := (3/(4 g p^4)) (g JJ[2 + q] - JJ[3 + q]) (g JJ[1 + q] - JJ[q]);
D202c[q_] := (5/(16 g p^6)) ((3 + 2 p^2) JJ[2 + q] - 6 g JJ[3 + q] + 3 JJ[4 + q]) ((3 + 2 p^2) JJ[1 + q] - 6 g JJ[q] + 3 JJ[q - 1]);
D020c[q_] := (5/(16 g p^6)) ((3 + 2 p^2) JJ[2 + q] - 6 g JJ[1 + q] + 3 JJ[q]) ((3 + 2 p^2) JJ[1 + q] - 6 g JJ[2 + q] + 3 JJ[3 + q]);
D121ref[q_] := (1/g) Fd[1, 2, q] Km[-q, 2, 1, bb];   (* = N16 D121, verified in harness 16 *)

(* ---------- PURE-PHYSICAL GENERATION (from D000 seed; no kernels for lp=0) ---------- *)
ClearAll[Da00, D00b];
Da00[0, q_] := D000c[q];  Da00[a_, q_] /; a < 0 := 0;
Da00[a_, q_] := Da00[a, q] = (1/Cm0[a]) ((g/p) Da00[a - 1, q] - (1/p) (JJ[1 + q]/JJ[2 + q]) Da00[a - 1, q + 1] - Cm0[a - 1] Da00[a - 2, q]);
D00b[0, q_] := D000c[q];  D00b[b_, q_] /; b < 0 := 0;
D00b[b_, q_] := D00b[b, q] = (1/Cm0[b]) ((g/p) D00b[b - 1, q] - (1/p) (JJ[2 + q]/JJ[1 + q]) D00b[b - 1, q - 1] - Cm0[b - 1] D00b[b - 2, q]);
Dh0l[h_, l_, q_] := Da00[h, q] D00b[l, q]/D000c[q];       (* P2 rank-one, lp=0 plane *)
(* middle climb (P3 step): D_{h,lp+1,0} from the two border ladders *)
mid[h_, lp_, q_] := Module[{Fh, F0},
   Fh = -(1/(Cm0[lp + 1] p)) (g Fd[h, lp, q] - Fd[h, lp, q - 1]) - (Cm0[lp]/Cm0[lp + 1]) Fd[h, lp - 1, q];
   F0 = (1/(Cm0[lp + 1] p)) (g Fd[lp, 0, q - 1] - Fd[lp, 0, q]) - (Cm0[lp]/Cm0[lp + 1]) Fd[lp - 1, 0, q - 1];
   ((-1)^(lp + 1)/g) Fh F0];
Dh20[h_, q_] := mid[h, 1, q];                            (* border D_{h,2,0} = climb lp 0->1->2 *)
D0lpl[l_, lp_, q_] := (* N1 reflection of the border *) If[lp == 2, Dh20[l, -3 - q], Dh0l[0, l, q]];
Dh2l[h_, l_, q_] := Dh20[h, q] D0lpl[l, 2, q]/Dh20[0, q]; (* P2 rank-one at lp=2 *)

qs = {1/2, -7/10};
report["C1  seed D000: pure-physical == closed form", Max@Table[Abs[Da00[0, q] - D000c[q]], {q, qs}]];
report["C2  D101 (lp=0): P0+rank-one == N8 closed form", Max@Table[Abs[Dh0l[1, 1, q] - D101c[q]], {q, qs}]];
report["C3  D202 (lp=0): P0+rank-one == N8 closed form", Max@Table[Abs[Dh0l[2, 2, q] - D202c[q]], {q, qs}]];
report["C4  D020 (lp=2): P3 middle climb == N8 closed form", Max@Table[Abs[Dh20[0, q] - D020c[q]], {q, qs}]];
report["C5  D121 (lp=2): P3 climb + N1 + rank-one == N16 closed form", Max@Table[Abs[Dh2l[1, 1, q] - D121ref[q]], {q, qs}]];

(* also: agreement of the two GENERATORS with each other (pure-physical vs the direct/lattice definition) *)
Ddir[l_, lp_, ldd_, q_] := (1/g) Kp[-1 - q, l, lp, bb] Km[-q, lp, ldd, bb];
report["C6  pure-physical == direct definition (spot: D101,D202,D020,D121)",
   Max[Abs[Dh0l[1, 1, 1/2] - Ddir[1, 0, 1, 1/2]], Abs[Dh0l[2, 2, 1/2] - Ddir[2, 0, 2, 1/2]],
       Abs[Dh20[0, 1/2] - Ddir[0, 2, 0, 1/2]], Abs[Dh2l[1, 1, 1/2] - Ddir[1, 2, 1, 1/2]]]];
Print["--- cross-check harness done ---"];
