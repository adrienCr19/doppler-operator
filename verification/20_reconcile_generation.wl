(* ::Package:: *)
(* 20_reconcile_generation.wl — reconcile the independent "Doppler Operator Generated from Its Seed"
   paper (doppler_generation_theorems.pdf) with this project's results.
   Notation: on nu^q, m=0,
     X(j,k; l,l1,l2,ldd) = (1/g) Kp[-1-j-q, l, l1] Km[-k-q, l2, ldd]     (off-diagonal family)
     D(j,k; l,lp,ldd)    = X(j,k; l,lp,lp,ldd)                            (physical, diagonal middle)
   Kp = kernel at +beta, Km at -beta, both m=0. *)

prec = 30;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{g = ga[b]}, 2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] /; l < 0 || lp < 0 := 0;
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
Cm0[l_] := If[l <= 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
JJ[w_, b_] := If[w == 0, 2 ArcTanh[b], With[{g = ga[b], p = ga[b] b}, ((g + p)^w - (g - p)^w)/w]];
report[name_, dev_, tol_: 10^-24] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (dev = ", ScientificForm[N[dev], 3], ")"];

Xo[j_, k_, q_, l_, l1_, l2_, ldd_, b_] := (1/ga[b]) Kp[-1 - j - q, l, l1, b] Km[-k - q, l2, ldd, b];
Dp[j_, k_, q_, l_, lp_, ldd_, b_] := Xo[j, k, q, l, lp, lp, ldd, b];

bb = 3/10; qv = 1/2;

(* ===== R1: gauge identity G2c:  (j+1,k+1)X(q) = (j,k)X(q+1) ===== *)
dev = Max@Table[Abs[Xo[j + 1, k + 1, qv, l, l1, l2, ldd, bb] - Xo[j, k, qv + 1, l, l1, l2, ldd, bb]],
   {j, {0, 1}}, {k, {-1, 0}}, {l, 0, 1}, {l1, 0, 1}, {l2, 0, 1}, {ldd, 0, 1}];
report["R1  gauge G2c: (j+1,k+1)X(q) = (j,k)X(q+1)", dev];

(* ===== R2: which rapidity ODE k-sign is correct? Ours (N14) vs paper (G5a) ===== *)
eta0 = ArcTanh[bb]; de = 10^-6;
fdD[j_, k_, l_, lp_, ldd_] := (Dp[j, k, qv, l, lp, ldd, Tanh[eta0 + de]] - Dp[j, k, qv, l, lp, ldd, Tanh[eta0 - de]])/(2 de);
(* our N14: Ohat -> -q, ldd terms carry -(ldd+k+Ohat)=-(ldd+k-q), +(ldd+1-k-Ohat)=+(ldd+1-k+q) *)
odeOurs[j_, k_, l_, lp_, ldd_] := -(l + 3 + j + qv) Cm0[l + 1] Dp[j, k, qv, l + 1, lp, ldd, bb] +
   (l - 2 - j - qv) Cm0[l] Dp[j, k, qv, l - 1, lp, ldd, bb] -
   (ldd + k - qv) Cm0[ldd + 1] Dp[j, k, qv, l, lp, ldd + 1, bb] +
   (ldd + 1 - k + qv) Cm0[ldd] Dp[j, k, qv, l, lp, ldd - 1, bb] - bb Dp[j, k, qv, l, lp, ldd, bb];
(* paper G5a: ldd terms carry -(ldd-k+Ohat)=-(ldd-k-q), +(ldd+1+k-Ohat)=+(ldd+1+k+q) *)
odePaper[j_, k_, l_, lp_, ldd_] := -(l + 3 + j + qv) Cm0[l + 1] Dp[j, k, qv, l + 1, lp, ldd, bb] +
   (l - 2 - j - qv) Cm0[l] Dp[j, k, qv, l - 1, lp, ldd, bb] -
   (ldd - k - qv) Cm0[ldd + 1] Dp[j, k, qv, l, lp, ldd + 1, bb] +
   (ldd + 1 + k + qv) Cm0[ldd] Dp[j, k, qv, l, lp, ldd - 1, bb] - bb Dp[j, k, qv, l, lp, ldd, bb];
devOurs = Max@Table[Abs[fdD[j, k, l, lp, ldd] - odeOurs[j, k, l, lp, ldd]], {j, {0, 1}}, {k, {1, -1}}, {l, 0, 1}, {lp, 0, 2, 2}, {ldd, 0, 1}];
devPaper = Max@Table[Abs[fdD[j, k, l, lp, ldd] - odePaper[j, k, l, lp, ldd]], {j, {0, 1}}, {k, {1, -1}}, {l, 0, 1}, {lp, 0, 2, 2}, {ldd, 0, 1}];
Print["R2  rapidity ODE k-sign at k = +-1 (FD-limited ~1e-8):"];
Print["     OUR N14 form:  dev = ", ScientificForm[N[devOurs], 3], "   ", If[TrueQ[N[devOurs] < 10^-7], "MATCHES FD", "FAILS"]];
Print["     PAPER G5a form: dev = ", ScientificForm[N[devPaper], 3], "   ", If[TrueQ[N[devPaper] < 10^-7], "MATCHES FD", "FAILS"]];

(* ===== R3: P1 outer pair relation (pure physical, companion eliminated) ===== *)
P1lhs[l_, lp_, ldd_] := Cm0[l + 1] Dp[0, 0, qv, l + 1, lp, ldd, bb] + Cm0[l] Dp[0, 0, qv, l - 1, lp, ldd, bb] -
   Cm0[ldd + 1] Dp[0, 0, qv + 1, l, lp, ldd + 1, bb] - Cm0[ldd] Dp[0, 0, qv + 1, l, lp, ldd - 1, bb];
P1rhs[l_, lp_, ldd_] := (ga[bb]/(ga[bb] bb)) (Dp[0, 0, qv, l, lp, ldd, bb] - Dp[0, 0, qv + 1, l, lp, ldd, bb]);
dev = Max@Table[Abs[P1lhs[l, lp, ldd] - P1rhs[l, lp, ldd]], {l, 0, 2}, {lp, 0, 2}, {ldd, 0, 2}];
report["R3  P1 outer pair relation (pure physical)", dev];

(* ===== R4: P0 axis ladder, explicit J-ratio coefficient ===== *)
P0lhs[a_] := Cm0[a + 1] Dp[0, 0, qv, a + 1, 0, 0, bb] + Cm0[a] Dp[0, 0, qv, a - 1, 0, 0, bb];
P0rhs[a_] := (ga[bb]/(ga[bb] bb)) Dp[0, 0, qv, a, 0, 0, bb] - (1/(ga[bb] bb)) (JJ[1 + qv, bb]/JJ[2 + qv, bb]) Dp[0, 0, qv + 1, a, 0, 0, bb];
dev = Max@Table[Abs[P0lhs[a] - P0rhs[a]], {a, 1, 3}];
report["R4  P0 axis ladder (D_a00), J-ratio coefficient", dev];
(* the a=1 seeded form *)
dev = Abs[Cm0[1] Dp[0, 0, qv, 1, 0, 0, bb] - ((ga[bb]/(ga[bb] bb)) Dp[0, 0, qv, 0, 0, 0, bb] - (1/(ga[bb] bb)) (JJ[1 + qv, bb]/JJ[2 + qv, bb]) Dp[0, 0, qv + 1, 0, 0, 0, bb])];
report["R4b P0 seeded: C1 D100(q) = (g/p)D000(q) - (1/p)(J1/J2)D000(q+1)", dev];

(* ===== R5: P2 rank-one ===== *)
dev = Max@Table[Abs[Dp[0, 0, qv, h1, lp, l1, bb] Dp[0, 0, qv, h2, lp, l2, bb] - Dp[0, 0, qv, h1, lp, l2, bb] Dp[0, 0, qv, h2, lp, l1, bb]],
   {lp, {0, 2}}, {h1, 0, 2}, {h2, 0, 2}, {l1, 0, 2}, {l2, 0, 2}];
report["R5  P2 rank-one: D_{h1 l' l1} D_{h2 l' l2} = D_{h1 l' l2} D_{h2 l' l1}", dev];

(* ===== R6: G5b mirror rapidity ODE (moves only the middle indices) ===== *)
fdX[j_, k_, l_, l1_, l2_, ldd_] := (Xo[j, k, qv, l, l1, l2, ldd, Tanh[eta0 + de]] - Xo[j, k, qv, l, l1, l2, ldd, Tanh[eta0 - de]])/(2 de);
(* paper G5b (Ohat -> -q): l1 terms +(l1-1-j+Ohat)=+(l1-1-j-q), -(l1+2+j-Ohat)=-(l1+2+j+q);
   l2 terms +(l2+2+k-Ohat)=+(l2+2+k+q), -(l2-1-k+Ohat)=-(l2-1-k-q) *)
mirror[j_, k_, l_, l1_, l2_, ldd_] := (l1 - 1 - j - qv) Cm0[l1 + 1] Xo[j, k, qv, l, l1 + 1, l2, ldd, bb] -
   (l1 + 2 + j + qv) Cm0[l1] Xo[j, k, qv, l, l1 - 1, l2, ldd, bb] +
   (l2 + 2 + k + qv) Cm0[l2 + 1] Xo[j, k, qv, l, l1, l2 + 1, ldd, bb] -
   (l2 - 1 - k - qv) Cm0[l2] Xo[j, k, qv, l, l1, l2 - 1, ldd, bb] - bb Xo[j, k, qv, l, l1, l2, ldd, bb];
dev = Max@Table[Abs[fdX[j, k, l, l1, l2, ldd] - mirror[j, k, l, l1, l2, ldd]], {j, {0, 1}}, {k, {0, 1}}, {l, 0, 1}, {l1, 0, 1}, {l2, 0, 1}, {ldd, 0, 1}];
report["R6  G5b mirror rapidity ODE (moves only middle indices)", dev, 10^-7];

(* ===== R7: F-bracket table spot checks vs our kernel columns =====
   F_{hk'}(q) = ^{-1-q}K_{hk'}(+b);  paper values in J-form: *)
Fpaper[0, 0] := JJ[2 + qv, bb]/(2 (ga[bb] bb));
Fpaper[1, 1] := With[{g = ga[bb], p = ga[bb] bb}, 3 (g JJ[1 + qv, bb] - (2 + p^2) JJ[2 + qv, bb] + g JJ[3 + qv, bb])/(2 p^3)];
Fpaper[2, 2] := With[{g = ga[bb], p = ga[bb] bb}, 5 (3 (3 + 2 p^2) JJ[qv, bb] - 12 g (3 + p^2) JJ[1 + qv, bb] + 2 (27 + 24 p^2 + 2 p^4) JJ[2 + qv, bb] - 12 g (3 + p^2) JJ[3 + qv, bb] + 3 (3 + 2 p^2) JJ[4 + qv, bb])/(8 p^5)];
Fours[h_, kp_] := Kp[-1 - qv, h, kp, bb];
dev = Max[Abs[Fpaper[0, 0] - Fours[0, 0]], Abs[Fpaper[1, 1] - Fours[1, 1]], Abs[Fpaper[2, 2] - Fours[2, 2]]];
report["R7  F-bracket table (F00, F11, F22) == our kernel columns ^{-1-q}K(+b)", dev];

(* ===== R8: bilinear assembly D_{hk'l} = ((-1)^{k'+l}/g) F_{hk'}(q) F_{k'l}(q-1), for D121 (new) ===== *)
F12q := Kp[-1 - qv, 1, 2, bb];
F21qm1 := Kp[-1 - (qv - 1), 2, 1, bb];
D121assembled := ((-1)^(2 + 1)/ga[bb]) F12q F21qm1;
D121direct := Dp[0, 0, qv, 1, 2, 1, bb];
report["R8  D121 bilinear assembly (paper G4b) == direct definition", Abs[D121assembled - D121direct]];

(* ===== R9: the 'new' elements the paper certifies (D212, D222, D303, D112): sanity (reflection N1a) ===== *)
diagset = {{2, 1}, {2, 2}, {3, 0}};
dev = Max@Table[Abs[Dp[0, 0, qv, ld[[1]], ld[[2]], ld[[1]], bb] - Dp[0, 0, -3 - qv, ld[[1]], ld[[2]], ld[[1]], bb]], {ld, diagset}];
report["R9  new diagonal elements (D212,D222,D303) reflection-invariant (N1a)", dev];

Print["--- reconciliation harness done ---"];
