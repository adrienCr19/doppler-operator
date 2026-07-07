(* ::Package:: *)
(* 08_paper_tables.wl — Check the user's rSZ paper (Mathematica/rSZ-boost/rSZ_paper.pdf) against
   the independently derived N8 closed forms:
   (P1) <D00> closed Bessel form [paper App. C]:
        <D00> = (K_{3-2O}(1/th) - K_1(1/th)) / (2 th (2-O)(1-O) K_2(1/th)),  O -> -q
   (P2) D101 channel data [paper App. I.0, sigma=2]: 5 terms (a_j, c_j), incl. corrected a_j=1 entry
   (P3) <D02> channel data [paper App. I.0, sigma=4]: 8 terms (a_j, c_j), D1 = (O-1)(O-2)
   Reconstruction claim: sinh^2(eta) cosh(eta) D(eta) = [Sum_j c_j cosh(a_j eta)]/sinh^sigma(eta),
   with Sum_j c_j = 0.  All D(eta) below are OUR closed forms (proofs/closed-forms.html). *)

Jw[w_, p_] := With[{eta = ArcSinh[p]}, If[w === 0, 2 eta, 2 Sinh[w eta]/w]];
D000c[q_, p_] := With[{g = Sqrt[1 + p^2]}, Jw[2 + q, p] Jw[1 + q, p]/(4 g p^2)];
D020c[q_, p_] := With[{g = Sqrt[1 + p^2]},
   5/(16 g p^6) ((3 + 2 p^2) Jw[2 + q, p] - 6 g Jw[1 + q, p] + 3 Jw[q, p]) *
                ((3 + 2 p^2) Jw[1 + q, p] - 6 g Jw[2 + q, p] + 3 Jw[3 + q, p])];
D101c[q_, p_] := With[{g = Sqrt[1 + p^2]},
   3/(4 g p^4) (g Jw[2 + q, p] - Jw[3 + q, p]) (g Jw[1 + q, p] - Jw[q, p])];
report[name_, dev_, tol_: 10^-20] := Print[name, ": ", If[dev < tol, "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

(* ---------- paper channel data (O = Ohat symbol) ---------- *)
tab101 = {
   {1,        3 (O^2 - 3 O + 6)/(8 O (O - 1) (O - 2) (O - 3))},
   {3,       -3/(8 O (O - 3))},
   {1 - 2 O, -3/(8 O (O - 1))},
   {3 - 2 O,  3 (O^2 - 3 O + 1)/(4 O (O - 1) (O - 2) (O - 3))},
   {5 - 2 O, -3/(8 (O - 2) (O - 3))}};
D1 = (O - 1) (O - 2);
tab02 = {
   {1,       -5 (O^4 - 6 O^3 + 17 O^2 - 24 O + 72)/(16 O (O - 1)^2 (O - 2)^2 (O - 3))},
   {3,        15 (O^2 - 3 O + 6)/(32 D1^2)},
   {5,       -5/(32 D1)},
   {1 - 2 O, -5 (2 O^3 - 12 O^2 + 13 O + 12)/(16 O D1^2)},
   {3 - 2 O,  15 (O^2 - 3 O - 3)/(16 D1^2)},
   {5 - 2 O, -5 (O + 1) (2 O^2 - 8 O + 3)/(16 (O - 3) D1^2)},
   {7 - 2 O,  5 O (O + 1)/(32 (O - 1) (O - 2)^2 (O - 3))},
   {1 + 2 O,  5 (O - 3) (O - 4)/(32 O (O - 1)^2 (O - 2))}};

(* ---------- P0: sum rules Sum_j c_j = 0 (symbolic) ---------- *)
Print["P0a  D101 table: Sum_j c_j = ", Simplify[Total[tab101[[All, 2]]]]];
Print["P0b  D02  table: Sum_j c_j = ", Simplify[Total[tab02[[All, 2]]]]];

(* ---------- P1: <D00> closed Bessel form vs thermal average of OUR D000 ---------- *)
fMB[p_?NumericQ, th_?NumericQ] := Exp[-Sqrt[1 + p^2]/th]/(th BesselK[2, 1/th]);
avg[Dfun_, q_?NumericQ, th_?NumericQ] := NIntegrate[p^2 fMB[p, th] Dfun[q, p], {p, 0, Infinity},
   WorkingPrecision -> 30, PrecisionGoal -> 14];
besselD00[q_, th_] := (BesselK[3 + 2 q, 1/th] - BesselK[1, 1/th])/(2 th (2 + q) (1 + q) BesselK[2, 1/th]);
dev = Max@Table[Abs[avg[D000c, q, th] - besselD00[q, th]],
   {q, {-7/10, 1/2, 13/10}}, {th, {1/50, 489/10000}}];  (* 10.2 keV and 25 keV *)
report["P1   <D000> == paper's closed Bessel form (App. C)", dev, 10^-12];

(* ---------- P2: D101 reconstruction identity at single momentum ---------- *)
(* claim: D101(eta) = [Sum_j c_j cosh(a_j eta)] / (sinh^4 eta cosh eta) *)
rec101[q_, eta_] := (Total[tab101 /. O -> -q /. {a_, c_} :> c Cosh[a eta]])/(Sinh[eta]^4 Cosh[eta]);
dev = Max@Table[Abs[D101c[q, Sinh[eta]] - rec101[q, eta]],
   {q, {-7/10, 1/2, 13/10}}, {eta, {3/10, 4/5, 3/2}}] // N[#, 20] &;
report["P2   D101(eta) == paper table reconstruction (sigma=2)", dev];

(* ---------- P3: D02 reconstruction identity at single momentum ---------- *)
rec02[q_, eta_] := (Total[tab02 /. O -> -q /. {a_, c_} :> c Cosh[a eta]])/(Sinh[eta]^6 Cosh[eta]);
dev = Max@Table[Abs[D020c[q, Sinh[eta]] - rec02[q, eta]],
   {q, {-7/10, 1/2, 13/10}}, {eta, {3/10, 4/5, 3/2}}] // N[#, 20] &;
report["P3   D020(eta) == paper table reconstruction (sigma=4)", dev];

(* ---------- P4: symbolic reconstruction (exact, all eta) ---------- *)
symdev[Dfun_, rec_, sig_] := Module[{lhs, rhs},
   lhs = Dfun[q, Sinh[eta]] Sinh[eta]^(sig + 2) Cosh[eta] // TrigToExp // Simplify;
   rhs = rec /. O -> -q /. {a_, c_} :> c Cosh[a eta] // Total // TrigToExp // Simplify;
   FullSimplify[lhs - rhs, Element[eta, Reals]]];
Print["P4a  symbolic: sinh^4 cosh * D101 - Sum c_j cosh(a_j eta) = ",
  symdev[D101c, tab101, 2]];
Print["P4b  symbolic: sinh^6 cosh * D020 - Sum c_j cosh(a_j eta) = ",
  symdev[D020c, tab02, 4]];

(* ---------- P5: thermally averaged channels via the paper's F-integrals vs direct averages ---------- *)
(* <D> = (z/K_2(z)) Integrate[e^{-z cosh eta} (Sum_j c_j cosh(a_j eta))/sinh^sigma eta] *)
avgF[tab_, sig_, q_?NumericQ, th_?NumericQ] := Module[{z = 1/th, integrand},
   integrand[eta_?NumericQ] = Exp[-z Cosh[eta]] (Total[tab /. O -> -q /. {a_, c_} :> c Cosh[a eta]])/Sinh[eta]^sig;
   (z/BesselK[2, z]) NIntegrate[integrand[eta], {eta, 0, Infinity}, WorkingPrecision -> 30, PrecisionGoal -> 14]];
dev = Max@Table[Abs[avg[D101c, q, th] - avgF[tab101, 2, q, th]], {q, {1/2, -7/10}}, {th, {489/10000}}];
report["P5a  <D101> direct == paper F-integral route", dev, 10^-12];
dev = Max@Table[Abs[avg[D020c, q, th] - avgF[tab02, 4, q, th]], {q, {1/2, -7/10}}, {th, {489/10000}}];
report["P5b  <D02>  direct == paper F-integral route", dev, 10^-12];

(* ---------- P6: S_th assembled from paper closed/Bessel pieces vs our exact average (07 harness anchor) ---------- *)
SthPaper[q_, th_] := besselD00[q, th] + avgF[tab02, 4, q, th]/10 - 1;
SthOurs[q_, th_] := avg[D000c, q, th] + avg[D020c, q, th]/10 - 1;
Print["P6   S_th(q=1/2, 25 keV): paper route = ", N[SthPaper[1/2, 489/10000], 14],
  "   our closed forms = ", N[SthOurs[1/2, 489/10000], 14]];
Print["--- paper-tables harness done ---"];
