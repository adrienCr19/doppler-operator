prec = 30;
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, m_, x_] := Sqrt[(2 l + 1)/(4 Pi) Factorial[l - m]/Factorial[l + m]] LegendreP[l, m, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, m_, b_] := Module[{g = ga[b]}, 2 Pi nint[Ybar[l, m, mu] Ybar[lp, m, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, m_, b_] /; l < Abs[m] || lp < Abs[m] := 0;
Km[d_?NumericQ, l_, lp_, m_, b_] := Km[d, l, lp, m, b] = KmRaw[d, l, lp, m, b];
e1 = 2/10; e2 = 3/10; Lmax = 16;
(* rapidity additivity (CR25 Eq. 19) for the kernel at -beta family: K(e1)K(e2) = K(e1+e2) *)
dev = Max@Table[Abs[Sum[Km[d, l, lp, m, Tanh[e1]] Km[d, lp, ldd, m, Tanh[e2]], {lp, Abs[m], Lmax}] - Km[d, l, ldd, m, Tanh[e1 + e2]]],
   {d, {0, 1}}, {m, 0, 1}, {l, m, 2}, {ldd, m, 2}];
Print["A1 rapidity additivity (CR25 19), Lmax=", Lmax, ": max dev = ", ScientificForm[N[dev], 3]];
