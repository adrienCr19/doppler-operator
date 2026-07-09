prec = 25;
sY[s_, l_, m_, mu_] := Module[{st2 = Sqrt[(1 - mu)/2], ct2 = Sqrt[(1 + mu)/2]},
   (-1)^(l + m - s) Sqrt[(2 l + 1)/(4 Pi) (l + m)! (l - m)!/((l + s)! (l - s)!)] *
    Sum[If[r + s - m < 0 || r > l - s || r + s - m > l + s, 0,
      Binomial[l - s, r] Binomial[l + s, r + s - m] (-1)^(l - r - s + m) *
       st2^(2 l - (2 r + s - m)) ct2^(2 r + s - m)], {r, Max[0, m - s], l - s}]];
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 14];
Do[Module[{mid = 2 Pi nint[mu sY[2, l, m, mu]^2]},
   Print["s=2 l=", l, " m=", m, ": <mu> = ", N[mid, 6], "   -sm/(l(l+1)) = ", N[-2 m/(l (l + 1)), 6],
     "   +sm/(l(l+1)) = ", N[2 m/(l (l + 1)), 6]]], {l, 2, 3}, {m, 0, 2}];
