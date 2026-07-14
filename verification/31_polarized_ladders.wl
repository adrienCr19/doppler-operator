(* ::Package:: *)
(* 31_polarized_ladders.wl — N24: the spin-weighted outer ladders (polarized N5/N6) and the
   ladder derivation of the polarized catalog.
   Mixed-spin lattice operator ((j,k) = weight lattice):
     (j,k)T(s1,s2)^m_{l lp ldd}(q) = (1/g) s1K^{-1-j-q,m}_{l lp}(+b) s2K^{-k-q,m}_{lp ldd}(-b)
   Ladders (sC^m_l = Sqrt[(l^2-m^2)(l^2-s^2)/(4l^2-1)]/l; sm-term signs sig5, sig6 pinned numerically):
     [N5s]  (j,k)T_l  = (1/(s1C^m_l p))[g (j,k)T_{l-1} - (j+1,k)T_{l-1}]
              + sig5 (s1 m/(l(l-1) s1C^m_l)) (j,k)T_{l-1} - (s1C^m_{l-1}/s1C^m_l) (j,k)T_{l-2}
     [N6s]  (j,k)T_ldd = (1/(s2C^m_ldd p))[g (j,k)T_{ldd-1} - (j,k-1)T_{ldd-1}]
              + sig6 (s2 m/(ldd(ldd-1) s2C^m_ldd)) (j,k)T_{ldd-1} - (s2C^m_{ldd-1}/s2C^m_ldd) (j,k)T_{ldd-2}
   Derivation: apply the verified spin-weighted kernel raising (P5 = RC26 B3 with our sm sign) to leg 1 at +b
   (row index = outgoing l; weight shift d -> d-1 is the lattice move j -> j+1), and via transposition
   (P3: s2K^{d}_{lp ldd}(-b) = s2K^{2-d}_{ldd lp}(+b)) to leg 2 (column index = incoming ldd; k -> k-1).
   Checks:
     Y1  the pinned ladders hold for spin pairs (2,2),(0,2),(2,0), m = 0,1,2, l/ldd = 3,4, two q, two b
     Y2  END-TO-END: seeds T^m_{222} (one integral each) -> N6s -> N5s -> T^m_{323} vs direct quadrature
     Y3  the ladder-climbed values == the N22 rationalization-lemma closed forms (route cross-check,
         the polarized analog of F11)  *)

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

(* lattice mixed operator; site convention: leg 1 at +beta = Ks[.., -b], leg 2 at -beta = Ks[.., +b] *)
TT[s1_, s2_, j_, k_, q_, l_, lp_, ldd_, m_, b_] :=
  (1/ga[b]) Ks[s1, -1 - j - q, l, lp, m, -b] Ks[s2, -k - q, lp, ldd, m, b];

(* ---------- ladders with sign parameters ---------- *)
N5s[sig_, s1_, s2_, j_, k_, q_, l_, lp_, ldd_, m_, b_] := With[{p = ga[b] b, C1 = Csm[s1, m, l], C0 = Csm[s1, m, l - 1]},
   (1/(C1 p)) (ga[b] TT[s1, s2, j, k, q, l - 1, lp, ldd, m, b] - TT[s1, s2, j + 1, k, q, l - 1, lp, ldd, m, b]) +
    sig (s1 m/(l (l - 1) C1)) TT[s1, s2, j, k, q, l - 1, lp, ldd, m, b] -
    (C0/C1) TT[s1, s2, j, k, q, l - 2, lp, ldd, m, b]];
N6s[sig_, s1_, s2_, j_, k_, q_, l_, lp_, ldd_, m_, b_] := With[{p = ga[b] b, C1 = Csm[s2, m, ldd], C0 = Csm[s2, m, ldd - 1]},
   (1/(C1 p)) (ga[b] TT[s1, s2, j, k, q, l, lp, ldd - 1, m, b] - TT[s1, s2, j, k - 1, q, l, lp, ldd - 1, m, b]) +
    sig (s2 m/(ldd (ldd - 1) C1)) TT[s1, s2, j, k, q, l, lp, ldd - 1, m, b] -
    (C0/C1) TT[s1, s2, j, k, q, l, lp, ldd - 2, m, b]];

(* pin the sm signs on one spin-2, m=1 configuration *)
pin[f_] := Module[{d1, d2},
  d1 = Abs[f[1] - TT[2, 2, 0, 0, 1/2, If[f === f5, 3, 2], 2, If[f === f5, 2, 3], 1, 3/10]];
  d2 = Abs[f[-1] - TT[2, 2, 0, 0, 1/2, If[f === f5, 3, 2], 2, If[f === f5, 2, 3], 1, 3/10]];
  If[d1 < d2, 1, -1]];
f5[sig_] := N5s[sig, 2, 2, 0, 0, 1/2, 3, 2, 2, 1, 3/10];
f6[sig_] := N6s[sig, 2, 2, 0, 0, 1/2, 2, 2, 3, 1, 3/10];
sig5 = pin[f5]; sig6 = pin[f6];
Print["Y0  sm-term signs pinned:  sig5 = ", sig5, ",  sig6 = ", sig6];

(* ---------- Y1: the ladders across the grid ---------- *)
dev = Max@Flatten@Table[Abs[N5s[sig5, ss[[1]], ss[[2]], 0, 0, q, l, 2, ldd, m, b] - TT[ss[[1]], ss[[2]], 0, 0, q, l, 2, ldd, m, b]],
   {ss, {{2, 2}, {0, 2}, {2, 0}}}, {q, {1/2, -7/10}}, {l, {3, 4}}, {ldd, {2, 3}}, {m, 0, 2}, {b, {3/10, 7/10}}];
report["Y1a [N5s] outgoing ladder: (2,2),(0,2),(2,0), m=0..2, l=3,4", dev];
dev = Max@Flatten@Table[Abs[N6s[sig6, ss[[1]], ss[[2]], 0, 0, q, l, 2, ldd, m, b] - TT[ss[[1]], ss[[2]], 0, 0, q, l, 2, ldd, m, b]],
   {ss, {{2, 2}, {0, 2}, {2, 0}}}, {q, {1/2, -7/10}}, {l, {2, 3}}, {ldd, {3, 4}}, {m, 0, 2}, {b, {3/10, 7/10}}];
report["Y1b [N6s] incoming ladder: (2,2),(0,2),(2,0), m=0..2, ldd=3,4", dev];
(* also on a shifted lattice point *)
dev = Max@Table[Abs[N5s[sig5, 2, 2, 0, -1, q, 3, 2, 2, m, 3/10] - TT[2, 2, 0, -1, q, 3, 2, 2, m, 3/10]], {q, {1/2}}, {m, 0, 2}];
report["Y1c [N5s] at lattice point (0,-1) (as used inside the climb)", dev];

(* ---------- Y2: end-to-end climb from the seeds ---------- *)
(* seeds: T^m_{222} on the lattice = products of the two seed kernel elements (one integral each);
   climb: N6s raises ldd 2->3 (needs (0,0) and (0,-1) seeds; the ldd-2=1 term vanishes for spin-2 legs),
          then N5s raises l 2->3 (needs (0,0) and (1,0) of the ldd=3 elements, themselves climbed). *)
climb323[s1_, s2_, q_, m_, b_] := Module[{T3 = Association[]},
  Do[T3[{j, k}] = N6s[sig6, s1, s2, j, k, q, 2, 2, 3, m, b], {j, 0, 1}, {k, 0, 0}];
  With[{p = ga[b] b, C1 = Csm[s1, m, 3], C0 = Csm[s1, m, 2]},
   (1/(C1 p)) (ga[b] T3[{0, 0}] - T3[{1, 0}]) + sig5 (s1 m/(3 2 C1)) T3[{0, 0}] -
    (C0/C1) TT[s1, s2, 0, 0, q, 1, 2, 3, m, b]]];
dev = Max@Flatten@Table[Abs[climb323[ss[[1]], ss[[2]], q, m, 3/10] - TT[ss[[1]], ss[[2]], 0, 0, q, 3, 2, 3, m, 3/10]],
   {ss, {{2, 2}, {0, 2}}}, {q, {1/2, -7/10}}, {m, 0, 2}];
report["Y2  end-to-end: seeds T^m_{222} -> N6s -> N5s -> T^m_{323} vs quadrature", dev];

(* ---------- Y3: ladder values == N22 closed forms (route cross-check) ---------- *)
mt = Null; (* J-reduction from harness 30 *)
Rv[s_, l_, m_, v_] := Module[{a = Abs[m + s], c = Abs[m - s], ex},
   ex = Simplify[sYv[s, l, m, Sqrt[vv], Sqrt[1 - vv]]/(vv^(a/2) (1 - vv)^(c/2)), Assumptions -> 0 < vv < 1];
   ex /. vv -> v];
JJ[w_, b_] := If[w === 0, 2 ArcTanh[b], With[{g = 1/Sqrt[1 - b^2], p = b/Sqrt[1 - b^2]}, ((g + p)^w - (g - p)^w)/w]];
Kclosed[s_, l_, lp_, m_, d_, b_] /; l < Max[Abs[m], Abs[s]] || lp < Max[Abs[m], Abs[s]] := 0;
Kclosed[s_, l_, lp_, m_, d_, b_] := Module[{a = Abs[m + s], c = Abs[m - s], K = Max[Abs[m], Abs[s]],
    v, vp, E1, P, kmin, cl, g = 1/Sqrt[1 - b^2]},
   v = (b + 1 - t)/(2 b); vp = v (1 - b)/t;
   E1 = Together[v^a (1 - v)^c Rv[s, l, m, v] Rv[s, lp, m, vp]];
   kmin = Exponent[Denominator[E1], t]; P = Expand[Numerator[E1] Cancel[t^kmin/Denominator[E1]]];
   cl = CoefficientList[P, t];
   2 Pi (1 - b)^(a/2) (1 + b)^(c/2)/b *
    Sum[cl[[jj + 1]] g^(-(jj - kmin - K) - 1) JJ[jj - kmin - K + 1 - d, b], {jj, 0, Length[cl] - 1}]];
TclosedN22[s1_, s2_, q_, l_, ldd_, m_, b_] := Sqrt[1 - b^2] Kclosed[s1, l, 2, m, -1 - q, -b] Kclosed[s2, 2, ldd, m, -q, b];
dev = Max@Flatten@Table[Abs[climb323[ss[[1]], ss[[2]], q, m, 3/10] - TclosedN22[ss[[1]], ss[[2]], q, 3, 3, m, 3/10]],
   {ss, {{2, 2}, {0, 2}}}, {q, {1/2}}, {m, 0, 2}];
report["Y3  ladder climb == N22 rationalization-lemma closed forms (T^m_{323})", dev];
Print["--- polarized ladders harness done ---"];
