(* ::Package:: *)
(* 13_generating.wl — N13: the l'-generating function of Doppler operators is the
   Henyey-Greenstein round trip (Funk-Hecke).
   Statement: for the round-trip operator with an axisymmetric rest-frame phase function
   p(cos Theta), the multipole matrix elements at fixed m are Sum_l' w_l' D^m_{l l' l''} where
   w_l' = the Funk-Hecke eigenvalue (Legendre coefficient) of p. For the HG phase function
     p_t(x) = (1/4pi) (1-t^2)/(1 + t^2 - 2 t x)^{3/2},   w_l' = t^l'   exactly.
   Position-space form used for the direct check (l = l'' = 0, m = 0; on nu^q):
     Sum_l' t^l' D_{0 l' 0}(q) =
       (1/(4 pi gamma)) Int dy1 dy2  p_t(y1.y2) D+(mu1)^(-3-q) D+(mu2)^(q)
   with D+(mu) = gamma (1 + beta mu), from the verified position-space representation of the
   two kernel factors (cf. N1 proof).  *)

prec = 20;
b = 3/10; g0 = 1/Sqrt[1 - b^2]; eta0 = ArcTanh[b]; qv = 1/2; tv = 2/5;
report[name_, dev_, tol_: 10^-8] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (dev = ", ScientificForm[N[dev], 3], ")"];

(* ---- G0: Funk-Hecke eigenvalues of the HG kernel: 2 pi Int p_t(x) P_l(x) dx = t^l ---- *)
pt[x_, t_] := (1/(4 Pi)) (1 - t^2)/(1 + t^2 - 2 t x)^(3/2);
dev = Max@Table[Abs[2 Pi NIntegrate[pt[x, tv] LegendreP[l, x], {x, -1, 1}, WorkingPrecision -> 25, PrecisionGoal -> 15] - tv^l], {l, 0, 6}];
report["G0  Funk-Hecke eigenvalues of HG kernel = t^l", dev, 10^-12];

(* ---- multipole side: Sum_l' t^l' D_{0 l' 0}(q) from kernel quadratures ---- *)
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> 30, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
KmRaw[d_?NumericQ, l_, lp_] := 2 Pi nint[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(g0 (1 + b mu))^d];
Km[d_?NumericQ, l_, lp_] := Km[d, l, lp] = KmRaw[d, l, lp];
Kp[d_?NumericQ, l_, lp_] := Kp[d, l, lp] = 2 Pi nint[Ybar[l, mu] Ybar[lp, (mu - b)/(1 - b mu)]/(g0 (1 - b mu))^d];
DD[lp_] := (1/g0) Kp[-1 - qv, 0, lp] Km[-qv, lp, 0];
Lmax = 16;
lhs = Sum[tv^lp DD[lp], {lp, 0, Lmax}];
Print["G1a multipole sum (Lmax=", Lmax, "): ", N[lhs, 12], "   tail est ~ ", ScientificForm[N[tv^(Lmax + 1) Abs[DD[Lmax]]], 2]];

(* ---- position-space side: 3D quadrature with the HG kernel ---- *)
Dp[mu_] := g0 (1 + b mu);
cosTh[m1_, m2_, dp_] := m1 m2 + Sqrt[(1 - m1^2) (1 - m2^2)] Cos[dp];
rhs = (2 Pi/(4 Pi g0)) NIntegrate[
    pt[cosTh[m1, m2, dp], tv] Dp[m1]^(-3 - qv) Dp[m2]^(qv),
    {m1, -1, 1}, {m2, -1, 1}, {dp, 0, 2 Pi},
    WorkingPrecision -> prec, PrecisionGoal -> 10, MaxRecursion -> 12];
Print["G1b HG round-trip quadrature: ", N[rhs, 12]];
report["G1  Sum_l' t^l' D_{0l'0} == HG round trip", Abs[lhs - rhs]/Abs[lhs], 10^-8];

(* ---- G2: Thomson consistency: the weights {1, 0, 1/10} reproduce the S_th combination ----
   (definitional check that the Funk-Hecke weights of the Thomson phase function
    p_T(x) = (3/16pi)(1+x^2) are w_0 = 1, w_1 = 0, w_2 = 1/10, w_{l>2} = 0) *)
pT[x_] := (3/(16 Pi)) (1 + x^2);
wl = Table[2 Pi NIntegrate[pT[x] LegendreP[l, x], {x, -1, 1}, WorkingPrecision -> 25, PrecisionGoal -> 15], {l, 0, 4}];
report["G2  Thomson Funk-Hecke weights = {1, 0, 1/10, 0, 0}", Max@Abs[wl - {1, 0, 1/10, 0, 0}], 10^-18];

Print["--- generating-function harness done ---"];
