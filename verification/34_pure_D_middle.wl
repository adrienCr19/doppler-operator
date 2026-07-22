(* ::Package:: *)
(* 34_pure_D_middle.wl — FULLY EXPLICIT pure-Doppler-operator recurrences for raising every index,
   with the middle index l' as the centerpiece. Every object appearing in every relation is a
   physical Doppler operator D_{l l' l''}(q) (possibly at a shifted spectral index q). NO kernels,
   NO bracket F, NO extended family X anywhere except in Ddir, which is the independent reference
   used only to check the pure-D relations. On nu^q, m=0. *)

prec = 40;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 20, MaxRecursion -> 24];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{gg = ga[b]},
   2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(gg (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] /; l < 0 || lp < 0 := 0;
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
CC[l_] := If[l <= 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
report[name_, dev_, tol_: 10^-14] := Print[name, ": ",
   If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (dev = ", ToString[N[dev], InputForm], ")"];

bb = 3/10; g = ga[bb]; p = ga[bb] bb;
JJ[w_] := If[w == 0, 2 ArcTanh[bb], ((g + p)^w - (g - p)^w)/w];

(* ---- independent reference (built from kernels): used ONLY to check ---- *)
Ddir[l_, lp_, ldd_, q_] := (1/g) Kp[-1 - q, l, lp, bb] Km[-q, lp, ldd, bb];
(* the seed closed form (allowed base case) *)
D000[q_] := JJ[2 + q] JJ[1 + q]/(4 g p^2);
report["SEED  D000 closed form vs direct", Abs[D000[1/2] - Ddir[0, 0, 0, 1/2]]];

(* =====================================================================
   OUTER INDEX (axis) — pure-D three-term ladders P0.  Raise l along the
   l''=l'=0 axis; and the mirror along l.  Every term a physical D.
   ===================================================================== *)
ClearAll[Da00, D00b];
Da00[0, q_] := D000[q];  Da00[a_, q_] /; a < 0 := 0;
Da00[a_, q_] := Da00[a, q] =
   (1/CC[a]) ((g/p) Da00[a - 1, q] - (1/p) (JJ[1 + q]/JJ[2 + q]) Da00[a - 1, q + 1] - CC[a - 1] Da00[a - 2, q]);
D00b[0, q_] := D000[q];  D00b[b_, q_] /; b < 0 := 0;
D00b[b_, q_] := D00b[b, q] =
   (1/CC[b]) ((g/p) D00b[b - 1, q] - (1/p) (JJ[2 + q]/JJ[1 + q]) D00b[b - 1, q - 1] - CC[b - 1] D00b[b - 2, q]);
report["P0a   axis ladder D_{a00} from D000 (a<=5, q in {1/2,-7/10})",
   Max@Table[Abs[Da00[a, q] - Ddir[a, 0, 0, q]], {a, 1, 5}, {q, {1/2, -7/10}}]];
report["P0b   axis ladder D_{00b} from D000 (b<=5, q in {1/2,-7/10})",
   Max@Table[Abs[D00b[b, q] - Ddir[0, 0, b, q]], {b, 1, 5}, {q, {1/2, -7/10}}]];

(* =====================================================================
   MIDDLE INDEX l'  —  THE CENTREPIECE.  Fully explicit pure-D border step:

   C_{lp+1} D_{h,lp+1,0}(q)
     = (D_{lp+1,00}(q-1)/p) [ g D_{h,lp,0}(q)/D_{lp,00}(q-1)
                              - (J_{q-1}/J_q) D_{h,lp,0}(q-1)/D_{lp-1... }... ]
   (see box below for the exact statement).  Inputs: border elements one and two
   middle-levels down (at q and q-1) and the axis elements D_{*,0,0} (from P0). *)
Dax[a_, q_] := If[a == 0, D000[q], Da00[a, q]];   (* axis element, pure-D via P0 *)

(* P3 middle-index raise, border l''=0, general outer h.  Pure D only. *)
midP3[h_, lp_, q_] :=
  (1/CC[lp + 1]) (
     (Dax[lp + 1, q - 1]/p) (
         g  DborderRef[h, lp, q]/Dax[lp, q - 1]
       - (JJ[q - 1]/JJ[q]) DborderRef[h, lp, q - 1]/Dax[lp, q - 2] )
     - If[lp >= 1, CC[lp] (Dax[lp + 1, q - 1]/Dax[lp - 1, q - 1]) DborderRef[h, lp - 1, q], 0] );

(* To test the RELATION in isolation, feed it the true lower border elements: *)
DborderRef[h_, lp_, q_] := Ddir[h, lp, 0, q];
report["P3rel middle border step vs direct (h<=3, lp=0..3, q in {1/2,-7/10})",
   Max@Table[Abs[midP3[h, lp, q] - Ddir[h, lp + 1, 0, q]],
             {h, 0, 3}, {lp, 0, 3}, {q, {1/2, -7/10}}]];

(* =====================================================================
   N1 reflection: supplies the OTHER border D_{0 l' l} from the first.
   D_{l l' l''}(q) = D_{l'' l' l}(-3 - q).
   ===================================================================== *)
report["N1    reflection D_{h l' 0}(q) = D_{0 l' h}(-3-q)",
   Max@Table[Abs[Ddir[h, lp, 0, q] - Ddir[0, lp, h, -3 - q]],
             {h, 0, 3}, {lp, 0, 3}, {q, {1/2, -7/10}}]];

(* =====================================================================
   P2 rank-one filling: any fixed-l' plane from its two borders.
   D_{h l' l}(q) = D_{h l' 0}(q) D_{0 l' l}(q) / D_{0 l' 0}(q).
   ===================================================================== *)
Dfill[h_, lp_, l_, q_] := Ddir[h, lp, 0, q] Ddir[0, lp, l, q]/Ddir[0, lp, 0, q];
report["P2    rank-one D_{h l' l} from borders (h,l<=3, lp=0..3)",
   Max@Table[Abs[Dfill[h, lp, l, 1/2] - Ddir[h, lp, l, 1/2]],
             {h, 0, 3}, {lp, 0, 3}, {l, 0, 3}]];

(* =====================================================================
   END-TO-END: generate a deep interior element D_{2,3,2} from the SEED D000
   using ONLY pure-D relations:  P0 (axes)  ->  P3 (raise middle on border, self
   feeding: each level built from the previous pure-D level)  ->  N1 (2nd border)
   ->  P2 (fill).  No kernels anywhere in the generation path.
   ===================================================================== *)
ClearAll[Dgen];               (* Dgen = the pure-D self-consistent border tower D_{h,l',0} *)
Dgen[h_, 0, q_] := Dax[h, q];                      (* middle 0  = axis (P0) *)
Dgen[h_, lp_, q_] /; lp < 0 := 0;
Dgen[h_, lp_, q_] := Dgen[h, lp, q] =
  (1/CC[lp]) (
     (Dax[lp, q - 1]/p) (
         g  Dgen[h, lp - 1, q]/Dax[lp - 1, q - 1]
       - (JJ[q - 1]/JJ[q]) Dgen[h, lp - 1, q - 1]/Dax[lp - 1, q - 2] )
     - If[lp >= 2, CC[lp - 1] (Dax[lp, q - 1]/Dax[lp - 2, q - 1]) Dgen[h, lp - 2, q], 0] );
(* border tower self-consistency (feeds its own output back in) *)
report["GENb  self-fed border tower D_{h l' 0} (h<=3, l'<=4)",
   Max@Table[Abs[Dgen[h, lp, 1/2] - Ddir[h, lp, 0, 1/2]], {h, 0, 3}, {lp, 0, 4}]];

(* second border via N1, then P2 fill — the full pure-D pipeline *)
DgenFill[h_, lp_, l_, q_] := Dgen[h, lp, q] Dgen[l, lp, -3 - q]/Dgen[0, lp, q];
report["GEND  END-TO-END D_{2,3,2}(1/2) from D000 (P0->P3->N1->P2, no kernels)",
   Abs[DgenFill[2, 3, 2, 1/2] - Ddir[2, 3, 2, 1/2]]];
report["GEND2 END-TO-END D_{3,2,1}(-7/10) from D000",
   Abs[DgenFill[3, 2, 1, -7/10] - Ddir[3, 2, 1, -7/10]]];

Print["--- pure-D middle-index harness done ---"];
