(* ::Package:: *)
(* 07_thermal_average.wl — Thermally average the NEW closed forms (N8) over the relativistic
   Maxwellian and compare against the established rSZ results:
     - moments <p^k>: CR26 (stag240) Eqs. 12-13
     - exact thermal SZ operator S_th = <D000 + D020/10> - 1: CR26 Eqs. 5-7, 11
     - asymptotic distortion: Y0, Y1 of Itoh, Kohyama & Nozawa (1998) (as also confirmed by CR26 Eq. 19)
   Everything below is built ONLY from the closed forms of proofs/closed-forms.html. *)

ga[b_] := 1/Sqrt[1 - b^2];
(* closed-form building blocks: J_w = 2 Sinh[w eta]/w, gamma = Cosh[eta], p = Sinh[eta] *)
Jw[w_, p_] := With[{eta = ArcSinh[p]}, If[w === 0, 2 eta, 2 Sinh[w eta]/w]];
D000c[q_, p_] := With[{g = Sqrt[1 + p^2]}, Jw[2 + q, p] Jw[1 + q, p]/(4 g p^2)];
D020c[q_, p_] := With[{g = Sqrt[1 + p^2]},
   5/(16 g p^6) ((3 + 2 p^2) Jw[2 + q, p] - 6 g Jw[1 + q, p] + 3 Jw[q, p]) *
                ((3 + 2 p^2) Jw[1 + q, p] - 6 g Jw[2 + q, p] + 3 Jw[3 + q, p])];
D101c[q_, p_] := With[{g = Sqrt[1 + p^2]},
   3/(4 g p^4) (g Jw[2 + q, p] - Jw[3 + q, p]) (g Jw[1 + q, p] - Jw[q, p])];
Sc[q_, p_] := D000c[q, p] + D020c[q, p]/10 - 1;

(* ---------- 1. symbolic p-expansion of the closed forms ---------- *)
serS = Collect[Normal@Series[Sc[q, p], {p, 0, 6}], p, Simplify];
f2 = Coefficient[serS, p, 2]; f4 = Coefficient[serS, p, 4]; f6 = Coefficient[serS, p, 6];
Print["T1  S(q,p) from closed forms: p^2: ", Simplify[f2], "   p^4: ", Simplify[f4]];
Print["    p^2 coeff == (q^2+3q)/3 (CR26 Eq.11): ", Simplify[f2 - (q^2 + 3 q)/3] === 0];
Print["    p^4 coeff == CR26 Eq.11 p^4 term:    ", Simplify[f4 - (42/5 (-q) + 7/2 q^2 - 21/5 (-q)^3 + 7/10 q^4)/15] === 0];

serD101 = Collect[Normal@Series[D101c[q, p], {p, 0, 4}], p, Simplify];
Print["T2  D101(q,p) closed form p-series: ", serD101];
Print["    matches ODE table (-(q+1)(q+2)/3 p^2 - (q^2+3q-4)(q^2+3q+2)/15 p^4): ",
  Simplify[serD101 - (-(q + 1) (q + 2)/3 p^2 - ((q^2 + 3 q - 4) (q^2 + 3 q + 2)/15) p^4)] === 0];

(* ---------- 2. relativistic Maxwellian moments (CR26 Eqs. 6, 12-13) ---------- *)
fMB[p_?NumericQ, th_?NumericQ] := Exp[-Sqrt[1 + p^2]/th]/(th BesselK[2, 1/th]);
mom[k_, th_] := mom[k, th] = NIntegrate[p^(k + 2) fMB[p, th], {p, 0, Infinity},
    WorkingPrecision -> 30, PrecisionGoal -> 15];
momCF[k_, th_] := 2 (2 th)^(k/2) BesselK[(k + 4)/2, 1/th] Gamma[(k + 3)/2]/(Sqrt[Pi] BesselK[2, 1/th]);
dev = Max@Table[Abs[mom[k, 1/50] - momCF[k, 1/50]], {k, {0, 2, 4, 6}}];
Print["T3  moment integrals vs CR26 Eq. 13 closed form: ", If[dev < 10^-20, "PASS", "FAIL"], " (", ScientificForm[N[dev], 3], ")"];
(* small-theta expansions used below: <p^2> = 3th(1 + 5th/2) + O(th^3), <p^4> = 15 th^2 + O(th^3) *)
r2[th_] := (mom[2, th] - (3 th + 15 th^2/2));
r4[th_] := (mom[4, th] - 15 th^2);
Print["T4  <p^2>-(3th+15th^2/2) ~ O(th^3): ratio at th=0.01/0.005 = ", N[r2[1/100]/r2[1/200], 4],
  " (expect ~8);   <p^4>-15th^2 ratio = ", N[r4[1/100]/r4[1/200], 4]];

(* ---------- 3. O(th^2) thermal average -> Itoh Y0, Y1 ---------- *)
(* <S>(q) = f2 <p^2> + f4 <p^4> + O(th^3) = th (3 f2) + th^2 (15 f4 + 15/2 f2) + O(th^3).
   Convert polynomial in q to the operator P(x d/dx) acting on n_pl and compare with Itoh. *)
Y0q = Expand[3 f2];
Y1q = Expand[15 f4 + 15/2 f2];
npl[x_] := 1/(Exp[x] - 1);
applyPoly[pol_, x_] := Module[{cl = CoefficientList[pol, q]},
   Sum[cl[[j + 1]] Nest[x D[#, x] &, npl[x], j], {j, 0, Length[cl] - 1}]];
Pf = x Exp[x]/(Exp[x] - 1)^2;
XT = x Coth[x/2]; ST = x/Sinh[x/2];
Y0itoh = XT - 4;
Y1itoh = -10 + (47/2) XT - (42/5) XT^2 + (7/10) XT^3 + ST^2 (-(21/5) + (7/5) XT);
Print["T5  Y0 from closed forms - Itoh Y0: ", FullSimplify[applyPoly[Y0q, x]/Pf - Y0itoh]];
Print["T6  Y1 from closed forms - Itoh Y1: ", FullSimplify[applyPoly[Y1q, x]/Pf - Y1itoh]];

(* ---------- 4. finite-temperature: exact <S>(q,th) from the closed forms ---------- *)
Savg[q_?NumericQ, th_?NumericQ] := NIntegrate[p^2 fMB[p, th] Sc[q, p], {p, 0, Infinity},
   WorkingPrecision -> 30, PrecisionGoal -> 13];
(* asymptotic prediction through th^2 and residual scaling ~ th^3 *)
asym[qv_, th_] := (f2 /. q -> qv) mom[2, th] + (f4 /. q -> qv) mom[4, th] + (f6 /. q -> qv) mom[6, th];
res[th_] := Abs[Savg[1/2, th] - asym[1/2, th]];
Print["T7  exact <S>(q=1/2) at th=0.02 (10.2 keV): ", N[Savg[1/2, 1/50], 12]];
Print["    asymptotic through <p^6>:              ", N[asym[1/2, 1/50], 12]];
Print["    residual(th=0.02)/residual(th=0.01) = ", N[res[1/50]/res[1/100], 4], "  (expect ~16 for O(p^8))"];

(* ---------- 5. thermally averaged D101 (kinematic dipole channel) ---------- *)
D101avg[q_?NumericQ, th_?NumericQ] := NIntegrate[p^2 fMB[p, th] D101c[q, p], {p, 0, Infinity},
   WorkingPrecision -> 30, PrecisionGoal -> 13];
serPred[q_, th_] := -(q + 1) (q + 2)/3 mom[2, th] - ((q^2 + 3 q - 4) (q^2 + 3 q + 2)/15) mom[4, th];
Print["T8  <D101>(q=1/2) exact at th=0.01: ", N[D101avg[1/2, 1/100], 12],
  "   vs p^4-series+moments: ", N[serPred[1/2, 1/100], 12]];
Print["--- thermal average harness done ---"];
