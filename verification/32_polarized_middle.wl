(* ::Package:: *)
(* 32_polarized_middle.wl — N24 (part 2): the polarized middle-index relation (M_s).
   Mixed-spin off-diagonal family ((j,k) = weight lattice; code convention: leg 1 = Ks[s1,.., -b],
   leg 2 = Ks[s2,.., +b], as in harness 31):
     (j,k)X(s1,s2)^m_{l; a,b; ldd}(q) = (1/g) s1K^{-1-j-q,m}_{l a}(+beta) s2K^{-k-q,m}_{b ldd}(-beta),
     T = X|_{a=b=lp}.
   Ingredients (three-term recursions for each middle slot, from the verified spin-weighted B10 = P5
   and transposition P3; sm signs e1, e2 pinned numerically):
     [A_s] leg-1 column at +beta:  s1K^{a}_{l,lp} = -(1/(C p))[g s1K^{a}_{l,lp-1} - s1K^{a+1}_{l,lp-1}]
             + e1 (s1 m/(lp(lp-1) C)) s1K^{a}_{l,lp-1} - (C'/C) s1K^{a}_{l,lp-2},   C = s1C^m_lp
     [B_s] leg-2 row at -beta: P5 as verified (weight companion lowered, k -> k+1 on the lattice).
   Solving each for the lp-element and multiplying (nine products; matched pairs reassemble):
     [M_s]  (j,k)T_{l,lp,ldd} = (1/(C1 C2 p^2)) [G1 G2 (j,k)T - G1 (j,k+1)T - G2 (j-1,k)T + (j-1,k+1)T]_{l,lp-1,ldd}
        + (C2'/(C1 C2 p)) [G1 (j,k)X - (j-1,k)X]_{l; lp-1,lp-2; ldd}
        + (C1'/(C1 C2 p)) [G2 (j,k)X - (j,k+1)X]_{l; lp-2,lp-1; ldd}
        + (C1' C2'/(C1 C2)) (j,k)T_{l,lp-2,ldd},
     with C1 = s1C^m_lp, C2 = s2C^m_lp, Ci' = siC^m_{lp-1}, Gi = g - ei si m p/(lp(lp-1)).
   Scalar limit s1 = s2 = 0: Gi -> g, C -> C^m: relation (M) of N10 recovered.
   Checks:
     V1  ingredient [A_s] (sign e1 pinned)      V2  ingredient [B_s] at the leg-2 argument (e2)
     V3  [M_s] vs quadrature: (2,2),(0,2),(2,0),(0,0), m = 0..2, lp = 3,4, two q, b up to 0.7
     V4  spin-2 collapse at lp -> 3: all lp-2 = 1 objects vanish, [M_s] closes on the T_{l,2,ldd} family
         alone — worked climb T(2,2)^m_{2,3,2} from the T^m_{222} lattice family  *)

prec = 25;
sYv[s_, l_, m_, st2_, ct2_] := (-1)^(m - s) Sqrt[(2 l + 1)/(4 Pi) (l + m)! (l - m)!/((l + s)! (l - s)!)] *
   Sum[If[r + s - m < 0 || r > l - s || r + s - m > l + s, 0,
     Binomial[l - s, r] Binomial[l + s, r + s - m] (-1)^(l - r - s + m) *
      st2^(2 l - (2 r + s - m)) ct2^(2 r + s - m)], {r, Max[0, m - s], l - s}];
sY[s_, l_, m_, mu_] := sYv[s, l, m, Sqrt[(1 - mu)/2], Sqrt[(1 + mu)/2]];
ga[b_] := 1/Sqrt[1 - b^2];
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 14, MaxRecursion -> 20];
KsRaw[s_, d_?NumericQ, l_, lp_, m_, b_] := Module[{g = ga[b]},
   2 Pi nint[sY[s, l, m, mu] sY[s, lp, m, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] /; l < Max[Abs[m], Abs[s]] || lp < Max[Abs[m], Abs[s]] := 0;
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] := Ks[s, d, l, lp, m, b] = KsRaw[s, d, l, lp, m, b];
Csm[s_, m_, l_] := If[l <= Max[Abs[m], Abs[s]] - 1 || l == 0, 0, Sqrt[(l^2 - m^2) (l^2 - s^2)/(4 l^2 - 1)]/l];
report[name_, dev_, tol_: 10^-18] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

XX[s1_, s2_, j_, k_, q_, l_, a_, b2_, ldd_, m_, b_] :=
  (1/ga[b]) Ks[s1, -1 - j - q, l, a, m, -b] Ks[s2, -k - q, b2, ldd, m, b];

(* ---------- V1: ingredient [A_s] — leg-1 column recursion (kernel at -b in code = +beta site) ---------- *)
Acand[e1_, s_, d_, l_, lp_, m_, b_] := With[{p = ga[b] b, C1 = Csm[s, m, lp], C0 = Csm[s, m, lp - 1]},
   -(1/(C1 p)) (ga[b] Ks[s, d, l, lp - 1, m, -b] - Ks[s, d + 1, l, lp - 1, m, -b]) +
    e1 (s m/(lp (lp - 1) C1)) Ks[s, d, l, lp - 1, m, -b] - (C0/C1) Ks[s, d, l, lp - 2, m, -b]];
e1 = If[Abs[Acand[1, 2, 1/2, 2, 3, 1, 3/10] - Ks[2, 1/2, 2, 3, 1, -3/10]] <
    Abs[Acand[-1, 2, 1/2, 2, 3, 1, 3/10] - Ks[2, 1/2, 2, 3, 1, -3/10]], 1, -1];
dev = Max@Flatten@Table[Abs[Acand[e1, s, d, l, lp, m, b] - Ks[s, d, l, lp, m, -b]],
   {s, {0, 2}}, {d, {1/2, -3/2}}, {l, {2, 3}}, {lp, {3, 4}}, {m, 0, 2}, {b, {3/10, 7/10}}];
report["V1  [A_s] column recursion at +beta (e1 = " <> ToString[e1] <> ")", dev];

(* ---------- V2: ingredient [B_s] — leg-2 row recursion (kernel at +b in code = -beta site) ---------- *)
Bcand[e2_, s_, d_, lp_, ldd_, m_, b_] := With[{p = ga[b] b, C1 = Csm[s, m, lp], C0 = Csm[s, m, lp - 1]},
   -(1/(C1 p)) (ga[b] Ks[s, d, lp - 1, ldd, m, b] - Ks[s, d - 1, lp - 1, ldd, m, b]) +
    e2 (s m/(lp (lp - 1) C1)) Ks[s, d, lp - 1, ldd, m, b] - (C0/C1) Ks[s, d, lp - 2, ldd, m, b]];
e2 = If[Abs[Bcand[1, 2, 1/2, 3, 2, 1, 3/10] - Ks[2, 1/2, 3, 2, 1, 3/10]] <
    Abs[Bcand[-1, 2, 1/2, 3, 2, 1, 3/10] - Ks[2, 1/2, 3, 2, 1, 3/10]], 1, -1];
dev = Max@Flatten@Table[Abs[Bcand[e2, s, d, lp, ldd, m, b] - Ks[s, d, lp, ldd, m, b]],
   {s, {0, 2}}, {d, {1/2, -3/2}}, {lp, {3, 4}}, {ldd, {2, 3}}, {m, 0, 2}, {b, {3/10, 7/10}}];
report["V2  [B_s] row recursion at -beta (e2 = " <> ToString[e2] <> ")", dev];

(* ---------- V3: the assembled [M_s] ---------- *)
Ms[s1_, s2_, j_, k_, q_, l_, lp_, ldd_, m_, b_] := Module[
  {p = ga[b] b, g = ga[b], C1, C2, C1p, C2p, G1, G2},
  C1 = Csm[s1, m, lp]; C2 = Csm[s2, m, lp]; C1p = Csm[s1, m, lp - 1]; C2p = Csm[s2, m, lp - 1];
  G1 = g - e1 s1 m p/(lp (lp - 1)); G2 = g - e2 s2 m p/(lp (lp - 1));
  (1/(C1 C2 p^2)) (G1 G2 XX[s1, s2, j, k, q, l, lp - 1, lp - 1, ldd, m, b]
     - G1 XX[s1, s2, j, k + 1, q, l, lp - 1, lp - 1, ldd, m, b]
     - G2 XX[s1, s2, j - 1, k, q, l, lp - 1, lp - 1, ldd, m, b]
     + XX[s1, s2, j - 1, k + 1, q, l, lp - 1, lp - 1, ldd, m, b]) +
   (C2p/(C1 C2 p)) (G1 XX[s1, s2, j, k, q, l, lp - 1, lp - 2, ldd, m, b]
     - XX[s1, s2, j - 1, k, q, l, lp - 1, lp - 2, ldd, m, b]) +
   (C1p/(C1 C2 p)) (G2 XX[s1, s2, j, k, q, l, lp - 2, lp - 1, ldd, m, b]
     - XX[s1, s2, j, k + 1, q, l, lp - 2, lp - 1, ldd, m, b]) +
   (C1p C2p/(C1 C2)) XX[s1, s2, j, k, q, l, lp - 2, lp - 2, ldd, m, b]];
dev = Max@Flatten@Table[Abs[Ms[ss[[1]], ss[[2]], 0, 0, q, l, lp, ldd, m, 3/10] -
     XX[ss[[1]], ss[[2]], 0, 0, q, l, lp, lp, ldd, m, 3/10]],
   {ss, {{2, 2}, {0, 2}, {2, 0}, {0, 0}}}, {q, {1/2, -7/10}}, {l, {2, 3}}, {lp, {3, 4}}, {ldd, {2}}, {m, 0, 2}];
report["V3a [M_s] vs quadrature: all spin pairs, m = 0..2, lp = 3,4 (b = 0.3)", dev];
dev = Max@Flatten@Table[Abs[Ms[2, 2, 0, 0, q, 2, 3, 2, m, 7/10] - XX[2, 2, 0, 0, q, 2, 3, 3, 2, m, 7/10]],
   {q, {1/2}}, {m, 0, 2}];
report["V3b [M_s] at b = 0.7 (spin-2 pair)", dev];

(* ---------- V4: the spin-2 collapse at lp = 3 ---------- *)
(* for s1 = s2 = 2 and lp = 3: every lp-2 = 1 object vanishes (kernel guard), so [M_s] closes on the
   lattice T_{l,2,ldd} family alone: the polarized quadrupole -> octupole step needs no off-diagonal X. *)
M3collapse[j_, k_, q_, l_, ldd_, m_, b_] := Module[
  {p = ga[b] b, g = ga[b], C1 = Csm[2, m, 3], C2 = Csm[2, m, 3], G1, G2},
  G1 = g - e1 2 m p/6; G2 = g - e2 2 m p/6;
  (1/(C1 C2 p^2)) (G1 G2 XX[2, 2, j, k, q, l, 2, 2, ldd, m, b]
     - G1 XX[2, 2, j, k + 1, q, l, 2, 2, ldd, m, b]
     - G2 XX[2, 2, j - 1, k, q, l, 2, 2, ldd, m, b]
     + XX[2, 2, j - 1, k + 1, q, l, 2, 2, ldd, m, b])];
dev = Max@Flatten@Table[Abs[M3collapse[0, 0, q, 2, 2, m, 3/10] - XX[2, 2, 0, 0, q, 2, 3, 3, 2, m, 3/10]],
   {q, {1/2, -7/10}}, {m, 0, 2}];
report["V4  spin-2 collapse: T(2,2)^m_{2,3,2} from the T^m_{2,2,2} lattice family alone", dev];
Print["--- polarized middle-index harness done ---"];
