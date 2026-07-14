(* ::Package:: *)
(* 28_dispersion.wl — N20: the Compton dispersion relation on the celestial line.
   Setting (N18/U5): on the spectral line q = -3/2 + i k the exact single-momentum Thomson matrix
     T(k)_{l l''} = D_{l0l''}(q) + D_{l2l''}(q)/10 - [delta - b mu-hat]_{l l''}
   is Hermitian with a persistent zero eigenvalue.  This harness:
     D1  IDENTIFIES the marginal mode: it is the boosted comoving equilibrium
         v_l(q) = 2pi Int Ybar_l(mu) [gamma(1 - b mu)]^q dmu   (the electron's rest-frame-isotropic
         power law, boosted to the lab) - per-row residuals at quadrature precision, two k, two beta;
         plus the wrong-exponent control (any other exponent fails at O(0.1)).
     D2  FLAT BAND: the zero eigenvalue persists at every k (the mode exists for every q on the line).
     D3  DISPERSION MAP: the damping branches lambda_j(k) at beta = 0.3 (l <= 8 truncation),
         connecting at beta -> 0 to the static Thomson rates (0, -1, -9/10, -1, -1, ...).
     D4  THERMAL AVERAGE: <T(k)> over the relativistic Maxwellian at 10 keV - Hermiticity survives,
         the zero mode does NOT (each p has its own comoving mode): the spectral gap g(k) = top
         eigenvalue < 0 is the exact per-mode thermalization rate; multiple scattering exponentiates,
         n_tau = e^{tau <T>} n_0.  *)

prec = 25;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 13, MaxRecursion -> 20];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{g = ga[b]}, 2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
Cm0[l_] := If[l == 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
report[name_, dev_, tol_: 10^-12] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

qof[k_] := -3/2 + I k;
Dop[q_, l_, lp_, ldd_, b_] := (1/ga[b]) Kp[-1 - q, l, lp, b] Km[-q, lp, ldd, b];
Tmat[k_, b_, LT_] := Table[Dop[qof[k], l, 0, ldd, b] + Dop[qof[k], l, 2, ldd, b]/10 -
    (KroneckerDelta[l, ldd] - b (Cm0[ldd + 1] KroneckerDelta[l, ldd + 1] + Cm0[ldd] KroneckerDelta[l, ldd - 1])),
   {l, 0, LT}, {ldd, 0, LT}];
vcom[k_, b_, LT_] := Table[2 Pi nintg[Ybar[l, mu] (ga[b] (1 - b mu))^qof[k]], {l, 0, LT}];

(* ---------- D1: identification of the marginal mode ---------- *)
(* rows l <= Lrow of the untruncated identity, with the source resolved to Lbig *)
Lrow = 5; Lbig = 14;
resrows[k_, b_] := Module[{v = vcom[k, b, Lbig]},
   Max@Table[Abs[Sum[(Dop[qof[k], l, 0, ldd, b] + Dop[qof[k], l, 2, ldd, b]/10 -
          (KroneckerDelta[l, ldd] - b (Cm0[ldd + 1] KroneckerDelta[l, ldd + 1] + Cm0[ldd] KroneckerDelta[l, ldd - 1]))) v[[ldd + 1]],
        {ldd, 0, Lbig}]]/Norm[v], {l, 0, Lrow}]];
report["D1a T(k).v = 0 for v = boosted comoving power law (k=0.7, beta=0.3)", resrows[7/10, 3/10], 10^-10];
report["D1b same at k = 1.3", resrows[13/10, 3/10], 10^-10];
report["D1c same at beta = 0.6 (k=0.7)", resrows[7/10, 6/10], 10^-8];
(* wrong-exponent control at k=0.7, beta=0.3, truncated matrix *)
TT8 = Tmat[7/10, 3/10, 8];
vex[e_] := Table[2 Pi nintg[Ybar[l, mu] (ga[3/10] (1 - (3/10) mu))^e], {l, 0, 8}];
ctrl = Min@Table[Norm[TT8 . vex[e]]/Norm[vex[e]], {e, {qof[7/10] + 1, qof[7/10] - 1, -3 - qof[7/10], Conjugate[qof[7/10]]}}];
Print["D1d control: best WRONG exponent gives |T.v|/|v| = ", ScientificForm[N[ctrl], 2], "  (identification is sharp)"];

(* ---------- D2: the flat band ---------- *)
flat = Max@Table[Min[Abs[Eigenvalues[N[Tmat[k, 3/10, 8]]]]], {k, {0, 1/2, 1, 2, 4}}];
report["D2  flat band: min |eigenvalue| = 0 at every k in {0, 0.5, 1, 2, 4}", flat, 10^-10];

(* ---------- D3: the dispersion map ---------- *)
Print["D3  Compton dispersion lambda_j(k) at beta = 0.3 (l <= 8; top 5 branches, descending):"];
Do[Module[{ev = Sort[Re /@ Eigenvalues[N[Tmat[k, 3/10, 8]]], Greater]},
   Print["    k = ", NumberForm[N[k], {4, 2}], ":  ", NumberForm[#, {7, 5}] & /@ Take[ev, 5]]],
  {k, {0, 1/2, 1, 2, 4}}];
ev0 = Sort[Re /@ Eigenvalues[N[Tmat[7/10, 1/100, 8]]], Greater];
Print["D3b beta -> 0 limit (beta = 0.01, k = 0.7): ", NumberForm[#, {7, 4}] & /@ Take[ev0, 5],
  "   [static Thomson rates 0, -9/10, -1, -1, ...]"];

(* ---------- D3c: the spectrum depends only on beta — scan at k = 0 ---------- *)
Print["D3c damping spectrum vs beta (k = 0, l <= 8; top 5):"];
Do[Module[{ev = Sort[Re /@ Eigenvalues[N[Tmat[0, b, 8]]], Greater]},
   Print["    beta = ", NumberForm[N[b], {4, 2}], ":  ", NumberForm[#, {7, 5}] & /@ Take[ev, 5]]],
  {b, {1/10, 3/10, 6/10}}];

(* ---------- D5: EXACT ISOSPECTRALITY — the unitary intertwiner W_k = [gamma(1-b mu)]^{ik} ----------
   Thomson scattering is elastic in the electron rest frame, so the full evolution commutes with
   multiplication by (nu')^{ik} = nu^{ik} [gamma(1-b mu)]^{ik}.  The nu^{ik} factor shifts the Mellin
   block k -> 0 while the angular factor W_k is unitary:  T(k) W_k = W_k T(0)  exactly.
   => the single-momentum dispersion relation is FLAT: spec T(k) = spec T(0) for all k
      (the ~1e-7 k-variation in D3 is pure l-truncation). *)
Wmat[k_, b_, L1_, L2_] := Table[2 Pi nintg[Ybar[l, mu] Ybar[ldd, mu] (ga[b] (1 - b mu))^(I k)], {l, 0, L1}, {ldd, 0, L2}];
Module[{k = 7/10, b = 3/10, Lr = 4, Lb = 14, WA, WB, TkW, WT0},
  WA = Wmat[k, b, Lb, Lr]; WB = Wmat[k, b, Lr, Lb];
  TkW = Table[Sum[(Dop[qof[k], l, 0, j, b] + Dop[qof[k], l, 2, j, b]/10 -
        (KroneckerDelta[l, j] - b (Cm0[j + 1] KroneckerDelta[l, j + 1] + Cm0[j] KroneckerDelta[l, j - 1]))) *
       WA[[j + 1, ldd + 1]], {j, 0, Lb}], {l, 0, Lr}, {ldd, 0, Lr}];
  WT0 = Table[Sum[WB[[l + 1, j + 1]] *
       (Dop[qof[0], j, 0, ldd, b] + Dop[qof[0], j, 2, ldd, b]/10 -
        (KroneckerDelta[j, ldd] - b (Cm0[ldd + 1] KroneckerDelta[j, ldd + 1] + Cm0[ldd] KroneckerDelta[j, ldd - 1]))),
      {j, 0, Lb}], {l, 0, Lr}, {ldd, 0, Lr}];
  report["D5  intertwining T(k) W_k = W_k T(0)  (k=0.7, beta=0.3, rows l<=4)", Max@Abs[TkW - WT0], 10^-8]];

(* ---------- D4: thermal average at 10 keV ---------- *)
th = 10/511;  (* theta_e *)
Needs["NumericalDifferentialEquationAnalysis`"];
gq = GaussianQuadratureWeights[14, 10^-4, 14/10, prec];
fMJ[p_] := p^2 Exp[-Sqrt[1 + p^2]/th];  (* un-normalized; normalized by the same quadrature *)
Znorm = Sum[gq[[i, 2]] fMJ[gq[[i, 1]]], {i, Length[gq]}];
LT4 = 6;
Tth[k_] := Sum[gq[[i, 2]] fMJ[gq[[i, 1]]] Tmat[k, gq[[i, 1]]/Sqrt[1 + gq[[i, 1]]^2], LT4], {i, Length[gq]}]/Znorm;
Do[Module[{M = Tth[k], evs},
   report["D4a <T(k)> Hermitian (k=" <> ToString[N[k]] <> ")", Max@Abs[M - ConjugateTranspose[M]], 10^-12];
   evs = Sort[Re /@ Eigenvalues[N[M]], Greater];
   Print["D4b <T(", NumberForm[N[k], {3, 1}], ")> spectrum (10 keV): top 4 = ",
     NumberForm[#, {9, 6}] & /@ Take[evs, 4]];
   Print["     spectral gap g(k) = ", ScientificForm[evs[[1]], 4],
     "  (thermalization rate per unit optical depth; n_tau = e^{tau<T>} n_0)"]],
  {k, {0, 7/10, 2, 4}}];
Print["--- dispersion harness done ---"];
