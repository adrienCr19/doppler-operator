(* ::Package:: *)
(* 27_eq20_structure.wl — structure of the RC26 Eq. 20 lemma (N15 open proof).
   Two conjectures that localize the mechanism:
     X0  EXCHANGE SYMMETRY of the spin-weighted harmonic theta-part:  sY_{lm} = mY_{ls}  (pointwise!)
         (Wigner little-d symmetry d^l_{-m,s} = d^l_{-s,m}; Goldberg phase (-1)^(m-s) is s<->m symmetric)
     X1  EXCHANGE SYMMETRY of the kernel:  d_sK^m_{llp}(b) = d_mK^s_{llp}(b)   (all d, l, lp, b)
         => any m-sum of kernel bilinears at spins (s1,s2) IS a spin-sum at azimuths (s1,s2):
            the "collective in m" mystery of Eq. 20 becomes a collective-spin statement.
     X2  BILINEAR SPIN-SWAP (same-argument form, general):
            W_{llp}(s1,s2;a,c) = Sum_m [ s1K^{a,m}_{llp}(+b) s2K^{c,m}_{llp}(+b) - (s1<->s2) ]  = 0 ?
         tested: (a) with the Eq.-20 weight constraint a+c=1;  (b) with INDEPENDENT weights;
                 (c) spin pairs (0,2),(0,1),(1,2);  (d) several (l,lp) incl. l != lp.
     X3  off-diagonal columns:  Sum_m s1K^{a,m}_{l lp1} s2K^{c,m}_{l lp2} - (s1<->s2),  lp1 != lp2
     X4  m-weighted diagnostic:  Sum_m m * [antisymmetrized product]  (source term of the ODE hierarchy)
     X5  original mixed-argument Eq.-20 shape but OFF the diagonal:  l != ldd, middle L = 2,3
   Note P10b (harness 15): per-m the identity FAILS, so whatever passes below is genuinely collective. *)

prec = 25;
sY[s_, l_, m_, mu_] := Module[{st2 = Sqrt[(1 - mu)/2], ct2 = Sqrt[(1 + mu)/2]},
   (-1)^(m - s) Sqrt[(2 l + 1)/(4 Pi) (l + m)! (l - m)!/((l + s)! (l - s)!)] *
    Sum[If[r + s - m < 0 || r > l - s || r + s - m > l + s, 0,
      Binomial[l - s, r] Binomial[l + s, r + s - m] (-1)^(l - r - s + m) *
       st2^(2 l - (2 r + s - m)) ct2^(2 r + s - m)], {r, Max[0, m - s], l - s}]];
ga[b_] := 1/Sqrt[1 - b^2];
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 14, MaxRecursion -> 20];
KsRaw[s_, d_?NumericQ, l_, lp_, m_, b_] := Module[{g = ga[b]},
   2 Pi nint[sY[s, l, m, mu] sY[s, lp, m, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] /; l < Max[Abs[m], Abs[s]] || lp < Max[Abs[m], Abs[s]] := 0;
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] := Ks[s, d, l, lp, m, b] = KsRaw[s, d, l, lp, m, b];
KsP[s_, d_?NumericQ, l_, lp_, m_, b_] := Ks[s, d, l, lp, m, -b];  (* kernel at +beta *)
report[name_, dev_, tol_: 10^-12] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

bb = 3/10;

(* ---------- X0: pointwise exchange of the harmonic theta-parts ---------- *)
dev = Max@Flatten@Table[Abs[sY[s, l, m, x] - sY[m, l, s, x]],
   {s, -2, 2}, {m, -2, 2}, {l, Max[2, 2], 4}, {x, {-9/10, -1/3, 1/10, 4/7, 24/25}}];
Print["X0  (diagnostic) UNSIGNED exchange sY_{lm} = mY_{ls}: deviates by ",
  ScientificForm[N[dev], 3], "  -> the sign (-1)^(s+m) in X6 is necessary; it cancels in kernels (X1)"];

(* ---------- X1: exchange symmetry of the kernel ---------- *)
dev = Max@Flatten@Table[Abs[KsP[s, d, l, lp, m, bb] - KsP[m, d, l, lp, s, bb]],
   {s, 0, 2}, {m, 0, 2}, {d, {7/10, -13/10}}, {l, 2, 3}, {lp, 2, 3}];
report["X1  kernel exchange  sK^{d,m} = mK^{d,s}  (s,m = 0..2)", dev];
dev = Max@Flatten@Table[Abs[KsP[s, d, l, lp, m, bb] - KsP[m, d, l, lp, s, bb]],
   {s, {-2, -1}}, {m, {0, 1, 2}}, {d, {7/10}}, {l, 2, 3}, {lp, 2, 3}];
report["X1b kernel exchange incl. negative s", dev];

(* ---------- X2: bilinear spin-swap, same-argument form ---------- *)
W[s1_, s2_, a_, c_, l_, lp_, b_] := Sum[
   KsP[s1, a, l, lp, m, b] KsP[s2, c, l, lp, m, b] - KsP[s2, a, l, lp, m, b] KsP[s1, c, l, lp, m, b],
   {m, -Min[l, lp], Min[l, lp]}];
(* (a) Eq.-20 constrained weights a + c = 1 (a = -1-q, c = 2+q), the transposed Eq.-20 shape *)
dev = Max@Flatten@Table[With[{a = -1 - q, c = 2 + q}, Abs[W[0, 2, a, c, l, lp, bb]]],
   {q, {1/2, -7/10}}, {l, 2, 4}, {lp, {2, 3}}];
report["X2a spin-swap, weights a+c=1 (Eq.-20 shape), (0,2), l,lp mixed", dev];
(* (b) INDEPENDENT weights *)
dev = Max@Flatten@Table[Abs[W[0, 2, 7/10, -13/10, l, lp, bb]], {l, 2, 4}, {lp, {2, 3}}];
report["X2b spin-swap, INDEPENDENT weights (a,c)=(0.7,-1.3), (0,2)", dev];
(* (c) other spin pairs *)
dev = Max@Flatten@Table[Abs[W[s1s2[[1]], s1s2[[2]], 7/10, -13/10, l, lp, bb]],
   {s1s2, {{0, 1}, {1, 2}}}, {l, 2, 3}, {lp, {2, 3}}];
report["X2c spin-swap, independent weights, pairs (0,1) and (1,2)", dev];
(* (d) different beta *)
dev = Abs[W[0, 2, 7/10, -13/10, 3, 2, 7/10]];
report["X2d spin-swap at beta = 0.7", dev];

(* ---------- X3: off-diagonal columns ---------- *)
W3[s1_, s2_, a_, c_, l_, lp1_, lp2_, b_] := Sum[
   KsP[s1, a, l, lp1, m, b] KsP[s2, c, l, lp2, m, b] - KsP[s2, a, l, lp1, m, b] KsP[s1, c, l, lp2, m, b],
   {m, -Min[l, lp1, lp2], Min[l, lp1, lp2]}];
devs = Table[Abs[W3[0, 2, 7/10, -13/10, l, 2, 3, bb]], {l, 2, 3}];
Print["X3  off-diagonal columns (lp1,lp2)=(2,3), l=2,3:  ", ScientificForm[N[devs], 3],
  "   (zero => hierarchy closes on antisymmetric objects)"];

(* ---------- X4: m-weighted diagnostic ---------- *)
M4[s1_, s2_, a_, c_, l_, lp_, b_] := Sum[
   m (KsP[s1, a, l, lp, m, b] KsP[s2, c, l, lp, m, b] - KsP[s2, a, l, lp, m, b] KsP[s1, c, l, lp, m, b]),
   {m, -Min[l, lp], Min[l, lp]}];
devs = Table[Abs[M4[0, 2, 7/10, -13/10, l, 2, bb]], {l, 2, 3}];
Print["X4  m-weighted antisym sum, l=2,3:  ", ScientificForm[N[devs], 3],
  "   (zero => ODE sources vanish by m-parity)"];

(* ---------- X5: original mixed-argument shape, off the outer diagonal ---------- *)
TT[s1_, s2_, q_, l_, lp_, ldd_, m_, b_] := (1/ga[b]) KsP[s1, -1 - q, l, lp, m, b] Ks[s2, -q, lp, ldd, m, b];
W5[q_, l_, lp_, ldd_, b_] := Sum[TT[0, 2, q, l, lp, ldd, m, b] - TT[2, 0, q, l, lp, ldd, m, b],
   {m, -Min[l, lp, ldd], Min[l, lp, ldd]}];
dev = Max@Flatten@Table[Abs[W5[q, l, lp, ldd, bb]], {q, {1/2}}, {l, {2, 3}}, {lp, {2, 3}}, {ldd, {2, 3}}];
report["X5  Eq.-20 shape off-diagonal: Sum_m T(0,2)_{l lp ldd} = Sum_m T(2,0)_{l lp ldd}", dev];

(* ---------- X6: sign-refined pointwise exchange:  sY_{lm} = (-1)^(s+m) mY_{ls} ---------- *)
dev = Max@Flatten@Table[Abs[sY[s, l, m, x] - (-1)^(s + m) sY[m, l, s, x]],
   {s, -2, 2}, {m, -2, 2}, {l, 2, 4}, {x, {-9/10, -1/3, 1/10, 4/7, 24/25}}];
report["X6  pointwise exchange with sign:  sY_{lm} = (-1)^(s+m) mY_{ls}", dev, 10^-15];

(* ---------- X7: ON-CURVE SWAP LEMMA — symbolic certificates ----------
   Half-angle rationalization: mu = (1-t^2)/(1+t^2), boost = dilation t -> la t
   (la = sqrt((1+b)/(1-b)); mu_- = (mu-b)/(1-b mu) <-> t_- = la t).
   Each sY(mu(t)) is RATIONAL in t, so the swap difference must cancel exactly (Together == 0).
   A certificate at given (l,lp,s1,s2) PROVES the bilinear spin-swap - and hence RC26 Eq. 20
   at that l - for ALL beta and ALL weights (a,c), to all orders. *)
mt[t_] := (1 - t^2)/(1 + t^2);
sYt[s_, l_, m_, t_] := Module[{st2 = t/Sqrt[1 + t^2], ct2 = 1/Sqrt[1 + t^2]},
   (-1)^(m - s) Sqrt[(2 l + 1)/(4 Pi) (l + m)! (l - m)!/((l + s)! (l - s)!)] *
    Sum[If[r + s - m < 0 || r > l - s || r + s - m > l + s, 0,
      Binomial[l - s, r] Binomial[l + s, r + s - m] (-1)^(l - r - s + m) *
       st2^(2 l - (2 r + s - m)) ct2^(2 r + s - m)], {r, Max[0, m - s], l - s}]];
Gt[s1_, s2_, l_, lp_, t_, u_] := Sum[
   sYt[s1, l, m, t] sYt[s1, lp, m, la t] sYt[s2, l, m, u] sYt[s2, lp, m, la u],
   {m, -Min[l, lp], Min[l, lp]}];
cert[s1_, s2_, l_, lp_] := Module[{d},
   d = Together[Gt[s1, s2, l, lp, t, u] - Gt[s1, s2, l, lp, u, t]];
   Print["X7  on-curve swap certificate (s1,s2)=(", s1, ",", s2, ") (l,lp)=(", l, ",", lp, "): ",
     If[d === 0, "PROVED (identically zero, all beta, all weights)", "NOT ZERO: " <> ToString[Short[d, 2]]]]];
Do[cert[0, 2, l, 2], {l, 2, 8}];
cert[0, 2, 3, 3]; cert[0, 2, 4, 3];
cert[0, 1, 5, 2]; cert[1, 2, 5, 2];

Print["--- Eq. 20 structure harness done ---"];
