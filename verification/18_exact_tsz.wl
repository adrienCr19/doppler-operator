(* ::Package:: *)
(* 18_exact_tsz.wl — N17: the exact (all-orders) thermal SZ operator, assembled.
     S_th(Ohat, th) = <D00> + <D02>/10 - 1
   with
     <D00>  = (K_{3+2q}(z) - K_1(z)) / (2 th (2+q)(1+q) K_2(z)),  z = 1/th   [Bessel-closed; q = -Ohat]
     <D02>  = (z/K_2(z)) Sum_j c_j(q) F(-2, a_j, z)                          [8-term channel data, verified
              symbolically against the N8 closed form in harness 08; NOT Bessel-closable (closure theorem)]
   Checks:
     X1  route A (Bessel + F-quadrature) == route B (Maxwell-Juttner average of the N8 closed forms)
         at 25 and 50 keV, q in {-4, -2, -0.7, 0.5, 1.3}   [two fully independent exact evaluations]
     X2  exact vs the p-Taylor tower (orders p^2..p^8 from the data product + exact Bessel moments):
         fractional tower errors tabulated -- the quantitative case for the exact operator
     X3  the F-part obeys its N9 characterization: cited (harness 09, L3b).  *)

prec = 30;
ga[b_] := 1/Sqrt[1 - b^2];
JJ[w_, b_] := If[w == 0, 2 ArcTanh[b], With[{g = ga[b], p = ga[b] b}, ((g + p)^w - (g - p)^w)/w]];
D000c[q_, b_] := With[{g = ga[b], p = ga[b] b}, JJ[2 + q, b] JJ[1 + q, b]/(4 g p^2)];
D020c[q_, b_] := With[{g = ga[b], p = ga[b] b},
   5/(16 g p^6) ((3 + 2 p^2) JJ[2 + q, b] - 6 g JJ[1 + q, b] + 3 JJ[q, b]) ((3 + 2 p^2) JJ[1 + q, b] - 6 g JJ[2 + q, b] + 3 JJ[3 + q, b])];
fMB[p_?NumericQ, th_?NumericQ] := Exp[-Sqrt[1 + p^2]/th]/(th BesselK[2, 1/N[th, 40]]);

(* channel data for <D02> (q-form; verified vs N8 closed forms, harness 08 P3/P4b) *)
tab02[q_] := With[{O = -q, D1 = (-q - 1) (-q - 2)}, {
    {1, -5 (O^4 - 6 O^3 + 17 O^2 - 24 O + 72)/(16 O (O - 1)^2 (O - 2)^2 (O - 3))},
    {3, 15 (O^2 - 3 O + 6)/(32 D1^2)},
    {5, -5/(32 D1)},
    {1 - 2 O, -5 (2 O^3 - 12 O^2 + 13 O + 12)/(16 O D1^2)},
    {3 - 2 O, 15 (O^2 - 3 O - 3)/(16 D1^2)},
    {5 - 2 O, -5 (O + 1) (2 O^2 - 8 O + 3)/(16 (O - 3) D1^2)},
    {7 - 2 O, 5 O (O + 1)/(32 (O - 1) (O - 2)^2 (O - 3))},
    {1 + 2 O, 5 (O - 3) (O - 4)/(32 O (O - 1)^2 (O - 2))}}];

(* route A: Bessel-closed monopole + one-quadrature F-part *)
besselD00[q_, th_] := With[{z = 1/N[th, 40]}, (BesselK[3 + 2 q, z] - BesselK[1, z])/(2 th (2 + q) (1 + q) BesselK[2, z])];
D02F[q_?NumericQ, th_?NumericQ] := With[{z = 1/N[th, 40]},
   (z/BesselK[2, z]) NIntegrate[Exp[-z Cosh[e]] (Total[tab02[q] /. {a_?NumericQ, c_} :> c Cosh[a e]])/Sinh[e]^4,
     {e, 0, Infinity}, WorkingPrecision -> prec, PrecisionGoal -> 14]];
SthA[q_?NumericQ, th_?NumericQ] := besselD00[q, th] + D02F[q, th]/10 - 1;

(* route B: direct Maxwell-Juttner average of the single-momentum closed forms *)
SthB[q_?NumericQ, th_?NumericQ] := NIntegrate[p^2 fMB[p, th] With[{b = p/Sqrt[1 + p^2]},
     D000c[q, b] + D020c[q, b]/10 - 1], {p, 0, Infinity}, WorkingPrecision -> prec, PrecisionGoal -> 13];

qsX = {-9/2, -4, -7/10, 1/2, 13/10}; ths = {489/10000, 978/10000};  (* 25 and 50 keV; q away from the removable poles q = 0,-1,-2,-3 *)
dev = Max@Table[Abs[SthA[q, th] - SthB[q, th]], {q, qsX}, {th, ths}];
Print["X1  exact S_th: Bessel+F route == Maxwell-Juttner route: ",
  If[TrueQ[N[dev] < 10^-11], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

(* ---- X2: exact vs the p-tower (orders 2..8) ---- *)
Cm0[l_] := If[l == 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
buildC[lp_] := Module[{c},
  c[k_][l_, ldd_] /; l < 0 || ldd < 0 := 0;
  c[0][l_, ldd_] := If[l == lp && ldd == lp, 1, 0];
  c[k_][l_, ldd_] := c[k][l, ldd] = Expand[(1/k) (
      -(l + 3 + q) Cm0[l + 1] c[k - 1][l + 1, ldd] + (l - 2 - q) Cm0[l] c[k - 1][l - 1, ldd]
      - (ldd - q) Cm0[ldd + 1] c[k - 1][l, ldd + 1] + (ldd + 1 + q) Cm0[ldd] c[k - 1][l, ldd - 1]
      - Sum[SeriesCoefficient[Tanh[x], {x, 0, j}] c[k - 1 - j][l, ldd], {j, 1, k - 1}])];
  c];
c0 = buildC[0]; c2 = buildC[2]; ord = 8;
toP[c_, l_, ldd_] := Collect[Normal@Series[Sum[c[k][l, ldd] et^k, {k, 0, ord}] /. et -> ArcSinh[pp], {pp, 0, ord}], pp, Expand];
Sser = toP[c0, 0, 0] + toP[c2, 0, 0]/10 - 1;
momN[k_, th_] := momN[k, th] = N[2 (2 th)^(k/2) BesselK[(k + 4)/2, 1/N[th, 40]] Gamma[(k + 3)/2]/(Sqrt[Pi] BesselK[2, 1/N[th, 40]]), 30];
tower[qq_, th_, N8ord_] := Module[{cl = PadRight[CoefficientList[Expand[Sser /. q -> qq], pp], ord + 1]},
   Sum[If[EvenQ[k] && k <= N8ord, cl[[k + 1]] momN[k, th], 0], {k, 2, N8ord}]];
Print["X2  exact vs tower (fractional error of the truncated tower):"];
Print["     th        q        exact S_th          p^4-tower err    p^8-tower err"];
Do[Module[{ex = SthA[q, th], t4, t8},
    t4 = tower[q, th, 4]; t8 = tower[q, th, 8];
    Print["    ", N[th, 3], "   ", N[q, 3], "   ", N[ex, 8],
      "   ", ScientificForm[N[Abs[(t4 - ex)/ex], 3], 2], "   ", ScientificForm[N[Abs[(t8 - ex)/ex], 3], 2]]],
  {th, ths}, {q, {-9/2, -4, 1/2}}];
Print["X3  (F-part ODE characterization: harness 09, check L3b — cited)"];
Print["--- exact tSZ harness done ---"];
