(* ::Package:: *)
(* 36_kernel_foundations.wl — the FOUNDATIONS under the Doppler-operator recurrences.
   Confirms, from the aberration-kernel integral definition, every ingredient the exact D-recurrences
   are built from, and then the reassembly of N5, N6 and the middle-index column recursion (A).
   Convention (matches the whole site / harness 34):
     Km[d,l,lp,b] = 2pi Int Ybar_l(mu) Ybar_lp(mu') [g(1+b mu)]^-d dmu,  mu'=(mu+b)/(1+b mu)
     site  D-kernel(-beta) = Km[.,.,.,+bb];  D-kernel(+beta) = Kp[.,.,.,bb] = Km[.,.,.,-bb].
   All checks are numeric quadrature vs the closed algebraic form claimed in the proof page. *)

prec = 40;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec,
   AccuracyGoal -> 26, PrecisionGoal -> 26, MaxRecursion -> 30];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{gg = ga[b]},
   2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(gg (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] /; l < 0 || lp < 0 := 0;
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
CC[l_] := If[l <= 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
report[name_, dev_, tol_: 10^-18] := Print[name, ": ",
   If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (dev = ", ToString[N[dev], InputForm], ")"];

bb = 3/10; g = ga[bb]; p = ga[bb] bb;
JJ[w_] := If[w == 0, 2 ArcTanh[bb], ((g + p)^w - (g - p)^w)/w];

(* ===== I0. the Legendre "mu" recurrence:  mu Ybar_l = C_{l+1} Ybar_{l+1} + C_l Ybar_{l-1} ===== *)
(* verified as 2pi Int mu Ybar_l Ybar_lp dmu = C_{l+1} delta_{lp,l+1} + C_l delta_{lp,l-1} *)
legMatrix[l_, lp_] := 2 Pi nintg[mu Ybar[l, mu] Ybar[lp, mu]];
report["I0  Legendre mu-matrix  <l|mu|lp> = C_{l+1}d(l+1)+C_l d(l-1)",
   Max@Table[Abs[legMatrix[l, lp] - (CC[l + 1] Boole[lp == l + 1] + CC[l] Boole[lp == l - 1])],
             {l, 0, 4}, {lp, 0, 5}]];

(* ===== I1. SEED:  Km[d,0,0,+-bb] = J_{1-d}/(2p)  ===== *)
report["I1  seed  Km[d,0,0,bb] = J_{1-d}/(2p)",
   Max@Table[Abs[Km[d, 0, 0, bb] - JJ[1 - d]/(2 p)], {d, {-1/2, 0, 3/5, -7/5}}]];
report["I1' seed sign-indep  Km[d,0,0,-bb] = J_{1-d}/(2p)",
   Max@Table[Abs[Km[d, 0, 0, -bb] - JJ[1 - d]/(2 p)], {d, {-1/2, 0, 3/5, -7/5}}]];

(* ===== I2. PARITY:  Km[d,l,lp,bb] = (-1)^{l+lp} Km[d,l,lp,-bb]  ===== *)
report["I2  parity  Km(bb) = (-1)^{l+lp} Km(-bb)",
   Max@Table[Abs[Km[d, l, lp, bb] - (-1)^(l + lp) Km[d, l, lp, -bb]],
             {d, {-1/2, 3/5}}, {l, 0, 3}, {lp, 0, 3}]];

(* ===== I3. TRANSPOSITION:  Km[d,l,lp,bb] = Km[2-d,lp,l,-bb]  (site: K_{l lp}(-b)=K^{2-d}_{lp l}(+b)) ===== *)
report["I3  transposition  Km[d,l,lp,bb] = Km[2-d,lp,l,-bb]",
   Max@Table[Abs[Km[d, l, lp, bb] - Km[2 - d, lp, l, -bb]],
             {d, {-1/2, 3/5, -7/5}}, {l, 0, 3}, {lp, 0, 3}]];

(* ===== I4. MASTER l-RAISING (KR), derived from Legendre + mu->weight, for Km (site -beta) =====
   C_l Km[d,l,lp,bb] = (1/p) Km[d-1,l-1,lp,bb] - (g/p) Km[d,l-1,lp,bb] - C_{l-1} Km[d,l-2,lp,bb]   *)
report["I4  master recursion (KR) for Km(-beta)",
   Max@Table[Abs[CC[l] Km[d, l, lp, bb]
      - ((1/p) Km[d - 1, l - 1, lp, bb] - (g/p) Km[d, l - 1, lp, bb] - CC[l - 1] Km[d, l - 2, lp, bb])],
      {d, {-1/2, 3/5, -7/5}}, {l, 1, 4}, {lp, 0, 3}]];
(* the +beta form (KR+) via parity: sign of the (l-1) bracket flips *)
report["I4+ (KR+) for Kp(+beta): D^d(+b) = +1/(C_l p)[g D^d_{l-1} - D^{d-1}_{l-1}] - (C_{l-1}/C_l) D^d_{l-2}",
   Max@Table[Abs[Kp[d, l, lp, bb]
      - ((1/(CC[l] p)) (g Kp[d, l - 1, lp, bb] - Kp[d - 1, l - 1, lp, bb]) - (CC[l - 1]/CC[l]) Kp[d, l - 2, lp, bb])],
      {d, {-1/2, 3/5, -7/5}}, {l, 1, 4}, {lp, 0, 3}]];

(* ===== the mu->weight identity used to prove (KR):  mu [g(1+b mu)]^-d = (1/b)[ g^-1 [g(1+bmu)]^-(d-1) - [g(1+bmu)]^-d ] ===== *)
report["I4a mu->weight identity (pointwise)",
   Max@Table[Module[{mm = N[m, 40], wd, wd1},
       wd = (g (1 + bb mm))^-d; wd1 = (g (1 + bb mm))^-(d - 1);
       Abs[mm wd - (1/bb) (g^-1 wd1 - wd)]], {d, {-1/2, 3/5}}, {m, {-7/10, 1/5, 9/10}}]];

(* ===== D-LEVEL REASSEMBLY: the exact recurrences, physical point (j,k)=(0,0) and lattice ===== *)
Dlat[j_, k_, l_, lp_, ldd_, q_] := (1/g) Kp[-1 - j - q, l, lp, bb] Km[-k - q, lp, ldd, bb];
Dp[l_, lp_, ldd_, q_] := Dlat[0, 0, l, lp, ldd, q];

(* N5 (raise l): D = 1/(C_l p)[g D_{l-1} - {}^{(1,0)}D_{l-1}] - (C_{l-1}/C_l) D_{l-2} *)
n5[j_, k_, l_, lp_, ldd_, q_] := ( (1/(CC[l] p)) (g Dlat[j, k, l - 1, lp, ldd, q] - Dlat[j + 1, k, l - 1, lp, ldd, q])
   - (CC[l - 1]/CC[l]) Dlat[j, k, l - 2, lp, ldd, q] );
report["R-N5  outgoing-index raising (lattice)",
   Max@Table[Abs[Dlat[j, k, l, lp, ldd, q] - n5[j, k, l, lp, ldd, q]],
      {q, {1/2, -7/10}}, {j, {0, 1}}, {k, {0, -1}}, {l, 1, 3}, {lp, 0, 2}, {ldd, 0, 2}]];

(* N6 (raise ldd): D = 1/(C_ldd p)[g D_{ldd-1} - {}^{(0,-1)}D_{ldd-1}] - (C_{ldd-1}/C_ldd) D_{ldd-2} *)
n6[j_, k_, l_, lp_, ldd_, q_] := ( (1/(CC[ldd] p)) (g Dlat[j, k, l, lp, ldd - 1, q] - Dlat[j, k - 1, l, lp, ldd - 1, q])
   - (CC[ldd - 1]/CC[ldd]) Dlat[j, k, l, lp, ldd - 2, q] );
report["R-N6  incoming-index raising (lattice)",
   Max@Table[Abs[Dlat[j, k, l, lp, ldd, q] - n6[j, k, l, lp, ldd, q]],
      {q, {1/2, -7/10}}, {j, {0, 1}}, {k, {0, -1}}, {l, 0, 2}, {lp, 0, 2}, {ldd, 1, 3}]];

(* Ingredient (A): COLUMN recursion of the +beta kernel (raise the shared middle index of factor 1)
   Kp[a,l,lp] = -1/(C_lp p)[g Kp[a,l,lp-1] - Kp[a+1,l,lp-1]] - (C_{lp-1}/C_lp) Kp[a,l,lp-2]   *)
report["R-A   column recursion of the +beta kernel (middle-index ingredient)",
   Max@Table[Abs[Kp[a, l, lp, bb]
      - (-(1/(CC[lp] p)) (g Kp[a, l, lp - 1, bb] - Kp[a + 1, l, lp - 1, bb]) - (CC[lp - 1]/CC[lp]) Kp[a, l, lp - 2, bb])],
      {a, {-3/2, -1/2, 4/5}}, {l, 0, 3}, {lp, 1, 3}]];

Print["--- kernel-foundations harness done ---"];
