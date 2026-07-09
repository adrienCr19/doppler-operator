(* ::Package:: *)
(* 19_sl2c.wl — N18: the Lorentz-group (SL(2,C)) dictionary.
   U1  boosts act on the celestial sphere as dilations of the stereographic coordinate:
       tan(theta'/2) = e^{-eta} tan(theta/2)  under  mu' = (mu+beta)/(1+beta mu)
   U2  the d = 1 aberration kernel is an ORTHOGONAL matrix (unitary point of the sphere sector):
       Sum_lp  1K_{l lp}(b) 1K_{l'' lp}(b) = delta_{l l''}      [row orthonormality]
   U3  general-d unitarity of the boost operator on the weight-matched photon measure:
       <X', Y'>_{nu^{1-2d} dnu dy} = <X, Y>  for weight-d fields boosted explicitly (2D quadrature)
   U4  the aberration kernel at COMPLEX Doppler weight Delta = 1 - i k (unitary principal series):
       (a) inversion:  Sum_lp K^{1-ik}_{l lp}(-b) K^{1-ik}_{lp l''}(+b) = delta
       (b) unitarity:  conj(K^{1-ik}_{lp l}(-b)) = [K^{1-ik}(-b)^{-1}]_{l lp} = K^{1-ik}_{l lp}(+b)
   U5  the exact Thomson single-momentum scattering matrix on the spectral line q = -3/2 + i k:
       T(k)_{l l''} = D_{l 0 l''} + D_{l 2 l''}/10 - [delta_{l l''} - b(C_{l''+1} d_{l,l''+1} + C_{l''} d_{l,l''-1})]
       is HERMITIAN (N1 at q-bar = -3-q); its (truncated) eigenvalues — the Compton "dispersion" — are real. *)

prec = 25;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 13, MaxRecursion -> 20];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{g = ga[b]},
   2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
Cm0[l_] := If[l == 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
report[name_, dev_, tol_: 10^-15] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

bb = 3/10; eta0 = ArcTanh[bb];

(* ---- U1: celestial dilation ---- *)
dev = Max@Table[With[{mup = (mu0 + bb)/(1 + bb mu0)},
    Abs[Sqrt[(1 - mup)/(1 + mup)] - Exp[-eta0] Sqrt[(1 - mu0)/(1 + mu0)]]], {mu0, {-9/10, -1/3, 0, 1/2, 95/100}}];
report["U1  boost = dilation: tan(th'/2) = e^-eta tan(th/2)", dev, 10^-20];

(* ---- U2: d=1 kernel is orthogonal (row orthonormality) ---- *)
Lmax = 16;
dev = Max@Table[Abs[Sum[Kp[1, l, lp, bb] Kp[1, ldd, lp, bb], {lp, 0, Lmax}] - KroneckerDelta[l, ldd]], {l, 0, 3}, {ldd, 0, 3}];
Print["U2  d=1 kernel orthogonality (Lmax=", Lmax, "): max dev = ", ScientificForm[N[dev], 3], " (truncation-limited)"];

(* ---- U3: invariance of the weight-matched inner product under an explicit boost ---- *)
(* weight-d fields: X'(nu',mu') = [1/(g(1+b mu'))]^d X(g nu'(1+b mu'), (mu'+b)/(1+b mu')) *)
(* test fields must vanish fast enough at nu -> 0 for the nu^(1-2d) measure (d up to 2) *)
XX[nu_, mu_] := nu^2 Exp[-nu] (1 + mu^2);  YY[nu_, mu_] := nu^3 Exp[-2 nu] (1 - mu/2);
pair[d_, F_, G_] := NIntegrate[F[nu, mu] G[nu, mu] nu^(1 - 2 d), {nu, 0, Infinity}, {mu, -1, 1},
   WorkingPrecision -> prec, PrecisionGoal -> 11];
boosted[d_, F_][nu_, mu_] := (1/(ga[bb] (1 + bb mu)))^d F[ga[bb] nu (1 + bb mu), (mu + bb)/(1 + bb mu)];
dev = Max@Table[Abs[pair[d, boosted[d, XX], boosted[d, YY]] - pair[d, XX, YY]], {d, {0, 1, 2}}];
report["U3  <X',Y'> = <X,Y> on nu^(1-2d) dnu dy (explicit boost, d = 0,1,2)", dev, 10^-9];

(* ---- U4: kernels at complex weight Delta = 1 - i k ---- *)
kk = 7/10; DD1 = 1 - I kk;
dev = Max@Table[Abs[Sum[Km[DD1, l, lp, bb] Kp[DD1, lp, ldd, bb], {lp, 0, Lmax}] - KroneckerDelta[l, ldd]], {l, 0, 3}, {ldd, 0, 3}];
Print["U4a inversion at Delta = 1 - ik (Lmax=", Lmax, "): max dev = ", ScientificForm[N[dev], 3], " (truncation-limited)"];
dev = Max@Table[Abs[Conjugate[Km[DD1, lp, l, bb]] - Kp[DD1, l, lp, bb]], {l, 0, 4}, {lp, 0, 4}];
report["U4b unitarity: conj(K^Delta(-b))^T = K^Delta(+b) at Delta = 1 - ik", dev, 10^-18];

(* ---- U5: the Thomson scattering matrix on the spectral line q = -3/2 + i k ---- *)
qk = -3/2 + I kk;
Dop[q_, l_, lp_, ldd_] := (1/ga[bb]) Kp[-1 - q, l, lp, bb] Km[-q, lp, ldd, bb];
Tmat[l_, ldd_] := Dop[qk, l, 0, ldd] + Dop[qk, l, 2, ldd]/10 -
   (KroneckerDelta[l, ldd] - bb (Cm0[ldd + 1] KroneckerDelta[l, ldd + 1] + Cm0[ldd] KroneckerDelta[l, ldd - 1]));
LT = 8; TT = Table[Tmat[l, ldd], {l, 0, LT}, {ldd, 0, LT}];
dev = Max@Abs[TT - ConjugateTranspose[TT]];
report["U5a T(k) Hermitian on the spectral line (per-element exact)", dev, 10^-18];
ev = Sort[Re /@ Eigenvalues[N[TT]], Greater];
Print["U5b Compton eigenvalues on the celestial line (k=0.7, beta=0.3, truncated at l<=", LT, "):"];
Print["     ", NumberForm[#, 6] & /@ Take[ev, 5]];
Print["     max |Im(eigenvalue)| = ", ScientificForm[Max[Abs[Im /@ Eigenvalues[N[TT]]]], 2], " (should be ~0 by Hermiticity)"];

(* U5c: does the marginal (zero) mode persist at another k? *)
kk2 = 13/10; qk2 = -3/2 + I kk2;
Tmat2[l_, ldd_] := (1/ga[bb]) Kp[-1 - qk2, l, 0, bb] Km[-qk2, 0, ldd, bb] +
   (1/(10 ga[bb])) Kp[-1 - qk2, l, 2, bb] Km[-qk2, 2, ldd, bb] -
   (KroneckerDelta[l, ldd] - bb (Cm0[ldd + 1] KroneckerDelta[l, ldd + 1] + Cm0[ldd] KroneckerDelta[l, ldd - 1]));
TT2 = Table[Tmat2[l, ldd], {l, 0, LT}, {ldd, 0, LT}];
ev2 = Sort[Re /@ Eigenvalues[N[TT2]], Greater];
Print["U5c top eigenvalues at k = 1.3: ", NumberForm[#, 6] & /@ Take[ev2, 3],
  "   (marginal mode persists: ", If[Abs[ev2[[1]]] < 10^-10, "YES", "NO"], ")"];
Print["--- SL(2,C) dictionary harness done ---"];
