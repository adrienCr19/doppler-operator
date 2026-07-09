(* ::Package:: *)
(* 15_polarized.wl — N15: the polarized (spin-weighted) sector.
   Spin-weighted harmonics via the Goldberg et al. (1967) formula (theta-part; azimuthal 2pi handled
   at matched m; theta-parts are real so kernels are real).
   Kernel (our convention, mirroring the verified s=0 one):
     K[s, d, l, lp, m](-b) = 2 pi Int sY_{lm}(mu') sY_{lp m}(mu+[mu']) / (g(1+b mu'))^d dmu'
   Checks:
     P1  orthonormality of the implemented sY
     P2  K(0) = delta
     P3  transposition (same s):  sK^d_{l lp}(-b) = sK^{2-d}_{lp l}(+b)
     P4  parity (s flips):        sK^d_{l lp}(-b) = (-1)^(l+lp) (-s)K^d_{l lp}(+b)
     P5  l-raising B10 with sm-term (RC26 Eq. B3) at s = 2, m = 0,1,2
     P6  seed integral (RC26 B5 first entry): 2K^0_22(-b) = (15/16) Int (1-mu^2)^2/[g(1+b mu)]^{d+2} dmu
     P7  polarized N1 (same spin):   sD(q)_{l lp ldd} = sD(-3-q)_{ldd lp l}
     P8  mixed-spin reflection:      T(s1,s2)(q)_{l lp ldd} = T(s2,s1)(-3-q)_{ldd lp l}
     P9  polarized parity:           T(s1,s2)(b) = (-1)^(l+ldd) T(-s1,-s2)(-b)
     P10 RC26 Eq. (20) identity, operator-exact: Sum_m T(0,2)^m_{l,2,l} = Sum_m T(2,0)^m_{l,2,l}
         for l = 2,3,4 at several q and beta (their check was a p-series; this is all-orders)
     P10b the same identity per single m (diagnostic: does it hold m-resolved?)
     P11 corollary: the m-summed mixed diagonals are invariant under q -> -3-q
     P12 kernel rapidity ODE with sm-term at s=2 (pin signs), column form  *)

prec = 25;
(* Goldberg formula, theta-part (e^{im phi} stripped); real *)
sY[s_, l_, m_, mu_] := Module[{st2 = Sqrt[(1 - mu)/2], ct2 = Sqrt[(1 + mu)/2]},
   (-1)^(m - s) Sqrt[(2 l + 1)/(4 Pi) (l + m)! (l - m)!/((l + s)! (l - s)!)] *
    Sum[If[r + s - m < 0 || r > l - s || r + s - m > l + s, 0,
      Binomial[l - s, r] Binomial[l + s, r + s - m] (-1)^(l - r - s + m) *
       st2^(2 l - (2 r + s - m)) ct2^(2 r + s - m)], {r, Max[0, m - s], l - s}]];
(* note: sin^{2l}(t/2) cot^{2r+s-m}(t/2) = sin^{2l-(2r+s-m)} cos^{2r+s-m}.
   Phase (-1)^(m-s): standard convention with <l|mu|l+1> = +sC^m_{l+1} and
   mu sY_lm = sC_{l+1} sY_{l+1,m} - (sm/(l(l+1))) sY_lm + sC_l sY_{l-1,m}  (measured). *)
ga[b_] := 1/Sqrt[1 - b^2];
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 14, MaxRecursion -> 20];
KsRaw[s_, d_?NumericQ, l_, lp_, m_, b_] := Module[{g = ga[b]},
   2 Pi nint[sY[s, l, m, mu] sY[s, lp, m, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] /; l < Max[Abs[m], Abs[s]] || lp < Max[Abs[m], Abs[s]] := 0;
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] := Ks[s, d, l, lp, m, b] = KsRaw[s, d, l, lp, m, b];
KsP[s_, d_?NumericQ, l_, lp_, m_, b_] := Ks[s, d, l, lp, m, -b];   (* kernel at +beta *)
Csm[s_, m_, l_] := If[l <= Max[Abs[m], Abs[s]] - 1 || l == 0, 0, Sqrt[(l^2 - m^2) (l^2 - s^2)/(4 l^2 - 1)]/l];
report[name_, dev_, tol_: 10^-18] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

bb = 3/10; qs = {1/2, -7/10};

(* P1 orthonormality *)
dev = Max@Flatten@Table[Abs[2 Pi nint[sY[2, l, m, mu] sY[2, lp, m, mu]] - KroneckerDelta[l, lp]], {m, 0, 2}, {l, 2, 4}, {lp, 2, 4}];
report["P1  orthonormality of implemented sY (s=2)", dev];

(* P2 identity at b=0 *)
dev = Max@Flatten@Table[Abs[Ks[2, d, l, lp, m, 0] - KroneckerDelta[l, lp]], {d, {0, 1}}, {m, 0, 2}, {l, 2, 4}, {lp, 2, 4}];
report["P2  sK(0) = delta", dev];

(* P3 transposition, same s *)
dev = Max@Flatten@Table[Abs[Ks[2, d, l, lp, m, bb] - KsP[2, 2 - d, lp, l, m, bb]], {d, {0, 1, -1/2}}, {m, 0, 2}, {l, 2, 4}, {lp, 2, 4}];
report["P3  sK^d_{llp}(-b) = sK^{2-d}_{lpl}(+b)", dev];

(* P4 parity with s-flip *)
dev = Max@Flatten@Table[Abs[Ks[2, d, l, lp, m, bb] - (-1)^(l + lp) KsP[-2, d, l, lp, m, bb]], {d, {0, 1}}, {m, 0, 2}, {l, 2, 4}, {lp, 2, 4}];
report["P4  sK^d(-b) = (-1)^(l+lp) (-s)K^d(+b)", dev];

(* P5 B10 with sm-term (RC26 B3):
   sK^d_{l lp}(-b) = -(1/(C p))[g sK^d_{l-1,lp} - sK^{d-1}_{l-1,lp}] - (sm/(l(l-1) C)) sK^d_{l-1,lp} - (C'/C) sK^d_{l-2,lp} *)
recB3[s_, d_, l_, lp_, m_, b_] := With[{g = ga[b], p = ga[b] b, C1 = Csm[s, m, l], C0 = Csm[s, m, l - 1]},
   -(1/(C1 p)) (g Ks[s, d, l - 1, lp, m, b] - Ks[s, d - 1, l - 1, lp, m, b]) +
    (s m/(l (l - 1) C1)) Ks[s, d, l - 1, lp, m, b] - (C0/C1) Ks[s, d, l - 2, lp, m, b]];
(* sm-term sign: + in OUR convention (kernel built from sY, <mu>_{sY} = -sm/(l(l+1)));
   RC26/CR25 print the opposite sign because their kernel uses (-s)Y (their footnote 4). *)
dev = Max@Flatten@Table[Abs[Ks[2, d, l, lp, m, bb] - recB3[2, d, l, lp, m, bb]], {d, {0, 1}}, {m, 0, 2}, {l, 3, 4}, {lp, 2, 3}];
report["P5  l-raising with sm-term (RC26 B3) at s=2", dev];

(* P6 seed integral (RC26 B5, m=0 entry) *)
dev = Max@Table[With[{g = ga[bb]},
    Abs[Ks[2, d, 2, 2, 0, bb] - (15/16) 2 Pi (1/(2 Pi)) nint[(1 - mu^2)^2/(g (1 + bb mu))^(d + 2)]]], {d, {0, 1, 2}}];
report["P6  2K^0_22(-b) = (15/16) Int (1-mu^2)^2 / D+^{d+2}", dev];

(* Doppler objects: same-spin and mixed-spin round trips on nu^q *)
TT[s1_, s2_, q_, l_, lp_, ldd_, m_, b_] := (1/ga[b]) KsP[s1, -1 - q, l, lp, m, b] Ks[s2, -q, lp, ldd, m, b];

(* P7 polarized N1, same spin *)
dev = Max@Flatten@Table[Abs[TT[2, 2, q, l, lp, ldd, m, bb] - TT[2, 2, -3 - q, ldd, lp, l, m, bb]],
   {q, qs}, {m, 0, 2}, {l, 2, 3}, {lp, 2, 3}, {ldd, 2, 3}];
report["P7  polarized N1: sD(q) = sD^T(-3-q)", dev];

(* P8 mixed-spin reflection *)
dev = Max@Flatten@Table[Abs[TT[0, 2, q, l, 2, ldd, m, bb] - TT[2, 0, -3 - q, ldd, 2, l, m, bb]],
   {q, qs}, {m, 0, 2}, {l, 2, 3}, {ldd, 2, 3}];
report["P8  mixed reflection: T(0,2)(q) = T(2,0)^T(-3-q)", dev];

(* P9 polarized parity *)
dev = Max@Flatten@Table[Abs[TT[2, 2, q, l, 2, ldd, m, bb] - (-1)^(l + ldd) TT[-2, -2, q, l, 2, ldd, m, -bb]],
   {q, {1/2}}, {m, 0, 2}, {l, 2, 3}, {ldd, 2, 3}];
report["P9  polarized parity: T(s1,s2)(b) = (-1)^(l+ldd) T(-s1,-s2)(-b)", dev];

(* P10 RC26 Eq. (20): m-summed mixed diagonals, (0,2) vs (2,0), operator-exact *)
msum[s1_, s2_, q_, l_, b_] := Sum[TT[s1, s2, q, l, 2, l, m, b], {m, -Min[l, 2], Min[l, 2]}];
Do[
  dev = Max@Table[Abs[msum[0, 2, q, l, b] - msum[2, 0, q, l, b]], {q, qs}, {b, {3/10, 7/10}}];
  report["P10 RC26 Eq. 20 at l=" <> ToString[l] <> ":  Sum_m T(0,2) = Sum_m T(2,0)", dev],
  {l, 2, 4}];

(* P10b per-m diagnostic *)
devs = Table[Abs[TT[0, 2, 1/2, 2, 2, 2, m, bb] - TT[2, 0, 1/2, 2, 2, 2, m, bb]], {m, 0, 2}];
Print["P10b per-m difference at l=2, m=0,1,2: ", ScientificForm[N[devs], 3], "  (identity may hold only m-summed)"];

(* P11 corollary: m-summed mixed diagonal invariant under q -> -3-q *)
dev = Max@Table[Abs[msum[0, 2, q, l, bb] - msum[0, 2, -3 - q, l, bb]], {q, qs}, {l, {2, 3}}];
report["P11 m-summed mixed diagonal invariant under q -> -3-q", dev];

(* P12 kernel ODE with sm-term at s=2 (column form at +beta), pin signs:
   d/deta sK^d_{l lp}(+b(eta)) = (lp+d) C_{lp+1} sK_{l,lp+1} - sm(1-d)/(lp(lp+1)) sK_{l,lp} - (lp+1-d) C_lp sK_{l,lp-1} *)
eta0 = ArcTanh[bb]; de = 10^-6;
fdK[s_, d_, l_, lp_, m_] := (KsP[s, d, l, lp, m, Tanh[eta0 + de]] - KsP[s, d, l, lp, m, Tanh[eta0 - de]])/(2 de);
ode[s_, d_, l_, lp_, m_] := (lp + d) Csm[s, m, lp + 1] KsP[s, d, l, lp + 1, m, bb] +
   s m (1 - d)/(lp (lp + 1)) KsP[s, d, l, lp, m, bb] - (lp + 1 - d) Csm[s, m, lp] KsP[s, d, l, lp - 1, m, bb];
(* same sm-sign flip as in recB3: our kernel carries sY, the papers' carries (-s)Y *)
dev = Max@Flatten@Table[Abs[fdK[2, d, l, lp, m] - ode[2, d, l, lp, m]], {d, {0, 1}}, {m, 0, 2}, {l, 2, 3}, {lp, 2, 3}];
report["P12 kernel ODE with sm-term at s=2 (column form, +b)", dev, 10^-8];

Print["--- polarized harness done ---"];
