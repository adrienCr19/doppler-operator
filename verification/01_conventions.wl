(* ::Package:: *)
(* 01_conventions.wl — Pin down aberration-kernel / boost-operator conventions
   and verify every literature identity numerically before deriving new ones.

   Conventions (s = 0 throughout this file):
     Kernel "at -beta"  :  Km[d, l, lp, m, b]  ==  ^dK^m_{l lp}(-beta)
        = Int dY' Ybar*_{lm}(mu') Ybar_{lp m}(mu[mu']) / (gamma(1+beta mu'))^d
        with mu[mu'] = (mu'+beta)/(1+beta mu')        [stag331 Eq. B4/B5]
     Kernel "at +beta"  :  Kp[d, l, lp, m, b] == ^dK^m_{l lp}(+beta) = Km with beta -> -beta.
     Boost operator on power laws nu^q:  Ohat nu^q = -q nu^q, so
        ^dBhat^m_{l lp}(nu, beta)  nu^q = Km[d - q, l, lp, m, beta] nu^q.
   Everything below prints PASS/FAIL with max abs deviations. *)

prec = 30;
nint[f_, opts___] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];

Ybar[l_, m_, x_] := Sqrt[(2 l + 1)/(4 Pi) Factorial[l - m]/Factorial[l + m]] LegendreP[l, m, x];

ga[b_] := 1/Sqrt[1 - b^2];

(* ^dK^m_{l lp}(-beta): azimuthal integral gives 2 Pi *)
KmRaw[d_?NumericQ, l_, lp_, m_, b_] := Module[{g = ga[b]},
   2 Pi nint[ Ybar[l, m, mu] Ybar[lp, m, (mu + b)/(1 + b mu)] / (g (1 + b mu))^d ]];
Km[d_?NumericQ, l_, lp_, m_, b_] /; l < Abs[m] || lp < Abs[m] := 0;
Km[d_?NumericQ, l_, lp_, m_, b_] := Km[d, l, lp, m, b] = KmRaw[d, l, lp, m, b];
Kp[d_?NumericQ, l_, lp_, m_, b_] := Km[d, l, lp, m, -b];

Cm[l_, m_] := If[l <= Abs[m], 0, Sqrt[(l^2 - m^2)/(4 l^2 - 1)]];

report[name_, dev_] := Print[name, ": ", If[dev < 10^-10, "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

bb = 3/10; pp = ga[bb] bb;

(* ---- Check 1: K(0) = identity ---- *)
dev = Max@Table[Abs[Km[d, l, lp, m, 0] - KroneckerDelta[l, lp]], {d, {0, 1, 2}}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}];
report["C1  K(beta=0) = delta", dev];

(* ---- Check 2: seed  ^dK^0_{00}(-b) = ((g+p)^(1-d) - (g-p)^(1-d)) / (2(1-d)p)   [stag240 Eq. 9a / B1] ---- *)
seed240[d_, b_] := With[{g = ga[b], p = ga[b] b}, ((g + p)^(1 - d) - (g - p)^(1 - d))/(2 (1 - d) p)];
dev = Max@Table[Abs[Km[d, 0, 0, 0, bb] - seed240[d, bb]], {d, {-2, -1, 1/2, 2, 3, 17/5}}];
report["C2  seed K00 closed form (stag240 B1)", dev];

(* ---- Check 3: seed  ^dK^m_{mm}(-b) = [(2m-1)!!]^2 (2m+1) / (2 (2m)!) * (1/g^{m+d}) 2F1((m+d)/2,(m+d+1)/2;3/2+m;p^2/g^2)
        [stag331 Eq. B4]  (note p^2/g^2 = b^2) ---- *)
seedB4[d_, m_, b_] := With[{g = ga[b]},
   ((2 m - 1)!!)^2 (2 m + 1)/(2 (2 m)!) 1/g^(m + d) Hypergeometric2F1[(m + d)/2, (m + d + 1)/2, 3/2 + m, b^2]];
dev = Max@Table[Abs[Km[d, m, m, m, bb] - seedB4[d, m, bb]], {d, {-1, 0, 1, 2, 7/3}}, {m, 0, 3}];
report["C3  seed Kmm 2F1 closed form (stag331 B4)", dev];

(* ---- Check 4: l-raising recursion [stag331 B3, s=0]:
   ^dK^m_{l lp}(-b) = -(1/C_l)[g ^dK_{l-1,lp} - ^{d-1}K_{l-1,lp}]/p - (C_{l-1}/C_l) ^dK_{l-2,lp}(-b) ---- *)
recB3[d_, l_, lp_, m_, b_] := With[{g = ga[b], p = ga[b] b},
   -(1/Cm[l, m]) (g Km[d, l - 1, lp, m, b] - Km[d - 1, l - 1, lp, m, b])/p - (Cm[l - 1, m]/Cm[l, m]) If[l - 2 >= Abs[m], Km[d, l - 2, lp, m, b], 0]];
dev = Max@Table[Abs[Km[d, l, lp, m, bb] - recB3[d, l, lp, m, bb]], {d, {0, 1, 2, -1}}, {m, 0, 1}, {l, Max[Abs[m], 1] + 1, 4}, {lp, Abs[m], 3}];
report["C4  l-raising recursion (stag331 B3 / stag240 B3)", dev];

(* ---- Check 5: transposition (stag698 13a):  ^dK^m_{lp l}(b) = ^{2-d}K^m_{l lp}(b) ---- *)
dev = Max@Table[Abs[Kp[d, lp, l, m, bb] - Kp[2 - d, l, lp, m, bb]], {d, {0, 1, 2, -1, 1/2}}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}];
report["C5  transposition K^d_{lp l}(b) = K^{2-d}_{l lp}(b) (stag698 13a)", dev];

(* ---- Check 6: parity (stag698 13b):  ^dK^m_{l lp}(-b) = (-1)^(l+lp) ^dK^m_{l lp}(b) ---- *)
dev = Max@Table[Abs[Km[d, l, lp, m, bb] - (-1)^(l + lp) Kp[d, l, lp, m, bb]], {d, {0, 1, 2, -1}}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}];
report["C6  parity (stag698 13b)", dev];

(* ---- Check 7: Doppler-weight raising (HLC26 A.7a, s=0):
   ^dK^m_{lp l}(b) = g ^{d+1}K^m_{lp l}(b) - g b [C_{lp+1} ^{d+1}K_{lp+1,l}(b) + C_{lp} ^{d+1}K_{lp-1,l}(b)] ---- *)
recA7a[d_, lp_, l_, m_, b_] := With[{g = ga[b]},
   g Kp[d + 1, lp, l, m, b] - g b (Cm[lp + 1, m] Kp[d + 1, lp + 1, l, m, b] + Cm[lp, m] If[lp - 1 >= Abs[m], Kp[d + 1, lp - 1, l, m, b], 0])];
dev = Max@Table[Abs[Kp[d, lp, l, m, bb] - recA7a[d, lp, l, m, bb]], {d, {0, 1, -1, 3/2}}, {m, 0, 1}, {lp, m, 3}, {l, m, 3}];
report["C7  Doppler-weight raising (HLC26 A.7a)", dev];

(* ---- Check 8: Doppler-weight lowering (HLC26 A.7b, s=0):
   ^dK^m_{lp l}(b) = g ^{d-1}K^m_{lp l}(b) + g b [C_{lp+1} ^{d-1}K_{lp+1,l}(b) + C_{lp} ^{d-1}K_{lp-1,l}(b)] ---- *)
recA7b[d_, lp_, l_, m_, b_] := With[{g = ga[b]},
   g Kp[d - 1, lp, l, m, b] + g b (Cm[lp + 1, m] Kp[d - 1, lp + 1, l, m, b] + Cm[lp, m] If[lp - 1 >= Abs[m], Kp[d - 1, lp - 1, l, m, b], 0])];
dev = Max@Table[Abs[Kp[d, lp, l, m, bb] - recA7b[d, lp, l, m, bb]], {d, {0, 1, -1, 3/2}}, {m, 0, 1}, {lp, m, 3}, {l, m, 3}];
report["C8  Doppler-weight lowering (HLC26 A.7b)", dev];

(* ---- Check 9: kernel ODE in rapidity (stag698 Eq. 28, s=0), column-shift form.
   Test BOTH candidate sign conventions against finite differences of Kp and Km. ---- *)
eta0 = ArcTanh[bb]; de = 10^-6;
fdKp[d_, l_, lp_, m_] := (Km[d, l, lp, m, -Tanh[eta0 + de]] - Km[d, l, lp, m, -Tanh[eta0 - de]])/(2 de);
fdKm[d_, l_, lp_, m_] := (Km[d, l, lp, m, Tanh[eta0 + de]] - Km[d, l, lp, m, Tanh[eta0 - de]])/(2 de);
ode28[K_, d_, l_, lp_, m_] := (lp + d) Cm[lp + 1, m] K[d, l, lp + 1, m, bb] - (lp + 1 - d) Cm[lp, m] If[lp - 1 >= Abs[m], K[d, l, lp - 1, m, bb], 0];
devP = Max@Table[Abs[fdKp[d, l, lp, m] - ode28[Kp, d, l, lp, m]], {d, {0, 1, 2}}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}];
devPneg = Max@Table[Abs[fdKp[d, l, lp, m] + ode28[Kp, d, l, lp, m]], {d, {0, 1, 2}}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}];
devM = Max@Table[Abs[fdKm[d, l, lp, m] - ode28[Km, d, l, lp, m]], {d, {0, 1, 2}}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}];
devMneg = Max@Table[Abs[fdKm[d, l, lp, m] + ode28[Km, d, l, lp, m]], {d, {0, 1, 2}}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}];
Print["C9  kernel rapidity ODE (stag698 28): FD tolerance ~1e-8"];
Print["     d/deta Kp = +ODE28[Kp]: ", ScientificForm[N[devP], 3], " | -ODE28[Kp]: ", ScientificForm[N[devPneg], 3]];
Print["     d/deta Km = +ODE28[Km]: ", ScientificForm[N[devM], 3], " | -ODE28[Km]: ", ScientificForm[N[devMneg], 3]];

(* ---- Check 10: inversion (stag698 16):  Sum_lp ^dK_{l lp}(-b) ^dK_{lp lpp}(b) = delta  (truncated sum) ---- *)
Lmax = 14;
dev = Max@Table[Abs[Sum[Km[d, l, lp, m, bb] Kp[d, lp, lpp, m, bb], {lp, Abs[m], Lmax}] - KroneckerDelta[l, lpp]],
    {d, {0, 1}}, {m, 0, 1}, {l, m, 2}, {lpp, m, 2}];
Print["C10 inversion (stag698 16), Lmax=", Lmax, ": max dev = ", ScientificForm[N[dev], 3], "  (limited by truncation)"];

(* ---- Check 11: stag240 Eq. (10a) closed form for Dhat_00 on power laws.
   Dhat_{0 lp 0} nu^q: Dhat = (1/g) * ^{-1}Bhat_{0lp}(nu,-b) ^0Bhat_{lp 0}(nu,b)
   On nu^q:  ^0Bhat_{lp 0}(nu, b) nu^q = Km[-q, lp, 0, 0, b] nu^q
             ^{-1}Bhat_{0 lp}(nu, -b) nu^q = Kp[-1-q, 0, lp, 0, b] nu^q
   Dhat00(q) = (1/g) Sum_lp Kp[-1-q,0,lp,0,b] Km[-q,lp,0,0,b]  ... but with lp fixed = 0 for the 000 element:
   D000(q) = (1/g) Kp[-1-q,0,0,0,b] Km[-q,0,0,0,b].
   stag240 (10a):  Dhat_00 = ((g+p)^{3-2Ohat} + (g-p)^{3-2Ohat} - 2g) / (4(2-Ohat)(1-Ohat) g p^2), Ohat -> -q *)
D000num[q_] := (1/ga[bb]) Kp[-1 - q, 0, 0, 0, bb] Km[-q, 0, 0, 0, bb];
D000cf[q_] := With[{g = ga[bb], p = pp}, ((g + p)^(3 + 2 q) + (g - p)^(3 + 2 q) - 2 g)/(4 (2 + q) (1 + q) g p^2)];
dev = Max@Table[Abs[D000num[q] - D000cf[q]], {q, {-7/10, 1/2, 13/10, 2}}];
report["C11 Dhat_000 closed form (stag240 10a)", dev];

Print["--- conventions harness done ---"];
