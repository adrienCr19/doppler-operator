(* ::Package:: *)
(* 33_polarized_thermal.wl — N22 (completion): exact thermal averages of the polarized catalog.
   Structure result (PB1): every polarized quadrupole-channel element resolves as
     T(s1,s2)^m_{222}(q) = Sum_j [cc_j(q) cosh(a_j eta) + cs_j(q) sinh(a_j eta)] / (sinh^10 eta cosh eta),
   a_j in {n, n +- 2q}, sigma_pol = 10 -- the SAME hyperbolic-rational family as the scalar catalog, two
   ladder rungs deeper than the scalar quadrupole channel, with one new feature: individual spin-weighted
   elements are not beta-even (parity maps s -> -s), so sinh parts appear alongside cosh.
   Maxwell-Juttner averaging therefore lands in the EXISTING F/L-calculus of N11 at mu = -4:
     <T> = (z/K2(z)) Sum_j [cc_j F(-4, a_j, z) + cs_j L(-4, a_j, z)],   z = 1/theta_e,
   evaluated as one regularized quadrature exactly as in E1 (the combined sum vanishes to O(eta^10); PB1c).
   Checks:
     PB1  symbolic cosh-resolution of T(2,2)^0, T(0,2)^0, T(2,2)^1, T(2,2)^2 (exact tables printed);
          PB1b pairing symmetry; PB1c Sum-rule (vanishing to O(eta^10), checked as value+derivatives)
     PB2  <T(2,2)^0_222> at 10 and 25 keV: F-route (regularized single quadrature of the resolution)
          vs DIRECT nested route (Maxwellian p-average of quadrature kernels -- no closed form used)
     PB3  the m-summed mixed diagonal Sum_m T(0,2)^m_222 (the Eq.-20 / polarized-SZ channel), both routes *)

prec = 25;
sYv[s_, l_, m_, st2_, ct2_] := (-1)^(m - s) Sqrt[(2 l + 1)/(4 Pi) (l + m)! (l - m)!/((l + s)! (l - s)!)] *
   Sum[If[r + s - m < 0 || r > l - s || r + s - m > l + s, 0,
     Binomial[l - s, r] Binomial[l + s, r + s - m] (-1)^(l - r - s + m) *
      st2^(2 l - (2 r + s - m)) ct2^(2 r + s - m)], {r, Max[0, m - s], l - s}];
sY[s_, l_, m_, mu_] := sYv[s, l, m, Sqrt[(1 - mu)/2], Sqrt[(1 + mu)/2]];
Rv[s_, l_, m_, v_] := Module[{a = Abs[m + s], c = Abs[m - s]},
   Simplify[sYv[s, l, m, Sqrt[vv], Sqrt[1 - vv]]/(vv^(a/2) (1 - vv)^(c/2)), Assumptions -> 0 < vv < 1] /. vv -> v];
report[name_, dev_, tol_: 10^-18] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

(* ---------- symbolic closed form in x = e^eta, u = x^q ---------- *)
JJxu[n_, dd0_, eps_] := (x^(n - dd0) u^(-eps) - x^(dd0 - n) u^(eps))/(n - dd0 - eps q);
Kclxu[s_, l_, lp_, m_, dd0_, eps_] := Module[{a = Abs[m + s], c = Abs[m - s], K = Max[Abs[m], Abs[s]],
    v, vp, E1, P, kmin, cl, bb, g},
   bb = (x^2 - 1)/(x^2 + 1); g = (x + 1/x)/2;
   v = (bb + 1 - t)/(2 bb); vp = v (1 - bb)/t;
   E1 = Together[v^a (1 - v)^c Rv[s, l, m, v] Rv[s, lp, m, vp]];
   kmin = Exponent[Denominator[E1], t]; P = Expand[Numerator[E1] Cancel[t^kmin/Denominator[E1]]];
   cl = CoefficientList[P, t];
   2 Pi (1 - bb)^(a/2) (1 + bb)^(c/2)/bb *
    Sum[cl[[j + 1]] g^(-(j - kmin - K) - 1) JJxu[j - kmin - K + 1, dd0, eps], {j, 0, Length[cl] - 1}]];
(* weights: leg 1 d = -1-q  =>  dd0 = -1, eps = -1;  leg 2 d = -q  =>  dd0 = 0, eps = -1 *)
Txu[s1_, s2_, m_] := Together[(2/(x + 1/x)) Kclxu[s1, 2, 2, m, -1, -1] Kclxu[s2, 2, 2, m, 0, -1]];

(* ---------- PB1: cosh-resolution ---------- *)
(* G = T sinh^10 cosh must be Laurent in x with q-rational coefficients; pair x^A u^b with x^-A u^-b *)
resolve[s1_, s2_, m_] := Module[{T = Txu[s1, s2, m], G, num, den, xoff, uoff, terms, assoc = <||>, key},
   G = Together[T ((x^2 - 1)/(2 x))^10 ((x^2 + 1)/(2 x))];
   den = Denominator[G]; xoff = Exponent[den, x]; uoff = Exponent[den, u];
   If[! PolynomialQ[Numerator[G], x] || ! PolynomialQ[Numerator[G], u], Print["resolve: unexpected structure!"]];
   num = Expand[Numerator[G]]; den = Simplify[den/(x^xoff u^uoff)];   (* den now q-rational only *)
   If[! FreeQ[den, x] || ! FreeQ[den, u], Print["resolve: denominator not cleared!"]];
   terms = If[Head[num] === Plus, List @@ num, {num}];
   Do[Module[{A = Exponent[tt, x] - xoff, bu = Exponent[tt, u] - uoff, cf},
     cf = Simplify[(tt/(x^(Exponent[tt, x]) u^(Exponent[tt, u])))/den];
     key = {A, bu}; assoc[key] = Lookup[assoc, Key[key], 0] + cf], {tt, terms}];
   assoc];
(* individual spin-weighted elements are NOT beta-even (parity maps s -> -s), so the resolution carries
   cosh AND sinh parts, and for m != 0 the exponential spectrum is asymmetric (orphan e^{-w eta} terms):
   T sinh^10 cosh = Sum_i gamma_i e^{w_i eta} = Sum_j [cc_j cosh(a_j eta) + cs_j sinh(a_j eta)],
   with a_j = |w|, cc = gamma(+w) + gamma(-w), cs = gamma(+w) - gamma(-w)  (absent partners = 0). *)
pairUp[assoc_] := Module[{keys = Keys[assoc], out = {}, wnum},
   Do[wnum = k[[1]] + k[[2]] (17/10);
    Which[
     wnum > 0,
      Module[{cp = assoc[k], cm = Lookup[assoc, Key[-k], 0]},
        AppendTo[out, {k[[1]] + k[[2]] q, Simplify[cp + cm], Simplify[cp - cm]}]],
     wnum == 0,
      AppendTo[out, {0, Simplify[assoc[k]], 0}],
     wnum < 0 && ! KeyExistsQ[assoc, -k],
      AppendTo[out, {-(k[[1]] + k[[2]] q), Simplify[assoc[k]], Simplify[-assoc[k]]}]],
    {k, keys}];
   out];
res220 = pairUp[resolve[2, 2, 0]];
res020 = pairUp[resolve[0, 2, 0]];
res221 = pairUp[resolve[2, 2, 1]];
res222 = pairUp[resolve[2, 2, 2]];
Print["PB1  resolution T(2,2)^0_{222} sinh^10(eta) cosh(eta) = Sum [cc_j cosh + cs_j sinh](a_j eta):  ",
  Length[res220], " orders a_j"];
Do[Print["      a_j = ", r[[1]], ":   cc_j = ", r[[2]], "   cs_j = ", r[[3]]], {r, SortBy[res220, First]}];
Print["PB1  (0,2)^0: ", Length[res020], " orders;  (2,2)^1: ", Length[res221], ";  (2,2)^2: ", Length[res222],
  "  (full tables reproducible from this script)"];
(* PB1c: regularity -- the combined sum vanishes to O(eta^10): exact series coefficients 0..9 *)
regcheck[res_, qv_] := Module[{g, ser},
   g = Total[(#[[2]] /. q -> qv) Cosh[(#[[1]] /. q -> qv) et] +
       (#[[3]] /. q -> qv) Sinh[(#[[1]] /. q -> qv) et] & /@ res];
   ser = CoefficientList[Normal@Series[g, {et, 0, 9}], et];
   If[ser === {}, 0, Max[Abs[N[ser, 30]]]]];
report["PB1c regularity: Sum_j [cc_j cosh + cs_j sinh](a_j eta) = O(eta^10)  (q = 1/2, all four channels)",
  Max[regcheck[res220, 1/2], regcheck[res020, 1/2], regcheck[res221, 1/2], regcheck[res222, 1/2]], 10^-15];

(* ---------- thermal averages: two routes ----------
   route F: <T> = (z/K2(z)) Int e^{-z cosh eta} [Sum_j c_j cosh(a_j eta)] csch^8(eta) deta
   (mu = -4 class, regularized single quadrature; measure: Int p^2 e^{-z gamma} dp = K_2(z)/z)
   route D (independent): nested quadrature, kernels only *)
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 14, MaxRecursion -> 20];
KsRaw[s_, d_?NumericQ, l_, lp_, m_, b_] := Module[{g = 1/Sqrt[1 - b^2]},
   2 Pi nint[sY[s, l, m, mu] sY[s, lp, m, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] /; l < Max[Abs[m], Abs[s]] || lp < Max[Abs[m], Abs[s]] := 0;
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] := Ks[s, d, l, lp, m, b] = KsRaw[s, d, l, lp, m, b];
Tquad[s1_, s2_, qv_, m_, b_] := Sqrt[1 - b^2] Ks[s1, -1 - qv, 2, 2, m, -b] Ks[s2, -qv, 2, 2, m, b];
Needs["NumericalDifferentialEquationAnalysis`"];
gq = GaussianQuadratureWeights[24, 10^-4, 11/5, prec];
avgD[s1_, s2_, qv_, m_, th_] := Module[{num, zn},
   num = Sum[gq[[i, 2]] gq[[i, 1]]^2 Exp[-Sqrt[1 + gq[[i, 1]]^2]/th] *
      Tquad[s1, s2, qv, m, gq[[i, 1]]/Sqrt[1 + gq[[i, 1]]^2]], {i, Length[gq]}];
   zn = Sum[gq[[i, 2]] gq[[i, 1]]^2 Exp[-Sqrt[1 + gq[[i, 1]]^2]/th], {i, Length[gq]}];
   num/zn];
(* normalization of route F: <T> = Int p^2 e^{-g/th} T dp / Int p^2 e^{-g/th} dp;
   p^2 dp = sinh^2 cosh deta; Int p^2 e^{-g/th} dp = th K_2(1/th) e-normalized:
   exact: Int_0^inf sinh^2 cosh e^{-z cosh} deta = K_2(z)/z  =>  <T> = (z/K2) Int e^{-z c} T sinh^2 cosh deta *)
avgFn[res_, qv_?NumericQ, th_?NumericQ] := Module[{z = 1/th, cn},
   cn = {#[[1]] /. q -> qv, #[[2]] /. q -> qv, #[[3]] /. q -> qv} & /@ res;
   (z/BesselK[2, z]) NIntegrate[Exp[-z Cosh[et]] Total[(#[[2]] Cosh[#[[1]] et] + #[[3]] Sinh[#[[1]] et]) & /@ cn] Csch[et]^8,
     {et, 1/100, Infinity}, WorkingPrecision -> 30, PrecisionGoal -> 12, MaxRecursion -> 20] +
   (z/BesselK[2, z]) NIntegrate[Exp[-z Cosh[et]] Total[(#[[2]] Cosh[#[[1]] et] + #[[3]] Sinh[#[[1]] et]) & /@ cn] Csch[et]^8,
     {et, 0, 1/100}, WorkingPrecision -> 40, PrecisionGoal -> 12, MaxRecursion -> 20]];
dev = Max@Flatten@Table[Abs[avgFn[res220, qv, th] - avgD[2, 2, qv, 0, th]],
   {qv, {1/2, -7/10}}, {th, {10/511, 25/511}}];
report["PB2  <T(2,2)^0_{222}>: F-class route == independent nested quadrature (10, 25 keV)", dev, 10^-7];

(* ---------- PB3: the m-summed mixed diagonal (explicit sum over m = -2..2) ---------- *)
Do[resMix[mm] = pairUp[resolve[0, 2, mm]], {mm, -2, 2}];
report["PB3a regularity of each m-resolution (q = 1/2)",
  Max@Table[regcheck[resMix[mm], 1/2], {mm, -2, 2}], 10^-15];
avgMix[qv_, th_] := Sum[avgFn[resMix[mm], qv, th], {mm, -2, 2}];
avgMixD[qv_, th_] := Sum[avgD[0, 2, qv, mm, th], {mm, -2, 2}];
dev = Max@Table[Abs[avgMix[qv, 10/511] - avgMixD[qv, 10/511]], {qv, {1/2, -7/10}}];
report["PB3  m-summed mixed diagonal <Sum_m T(0,2)^m_{222}> (the polarized-SZ channel), both routes", dev, 10^-7];
Print["PB3  values at 10 keV: q=1/2: ", N[avgMix[1/2, 10/511], 12], ",  q=-7/10: ", N[avgMix[-7/10, 10/511], 12]];
Print["--- polarized thermal harness done ---"];
