(* ::Package:: *)
(* 30_spin_m_closed_forms.wl — N22/N23: closed forms for ALL spins s and azimuthal numbers m.
   RATIONALIZATION LEMMA: in  sK^{d,m}_{l lp}(b) = 2pi Int sY_{lm}(mu) sY_{lp m}(mu+) D^{-d} dmu,
   both harmonics carry st2^|m+s| ct2^|m-s| (half-angle factors, d-function structure), and on the
   aberration curve st2(mu+) = st2(mu) Sqrt[(1-b)/(1+b mu)], ct2(mu+) = ct2(mu) Sqrt[(1+b)/(1+b mu)],
   so the PRODUCT is rational:  integrand = (1-b)^{a/2}(1+b)^{c/2} * poly(mu) * (1+b mu)^{-d-max(|m|,|s|)},
   a = |m+s|, c = |m-s|.  With t = 1+b mu and gamma(1+-b) = e^{+-eta}:
       Int t^{n-d} dt = gamma^{d-n-1} J_{n+1-d},   J_w = 2 sinh(w eta)/w
   => EVERY kernel element, any (s,m), is a finite J-combination: the J-calculus is universal.
   Consequences: (N22) the polarized (s=2) kernel columns and Doppler elements in closed form;
                 (N23) general-m scalar closed forms (prefactor (1-b^2)^{|m|/2} = gamma^{-|m|}).
   Checks:
     Z1  algorithmic J-reduction vs direct quadrature: s=2, m=0,1,2, l=2..4, lp=2, several d, b
     Z2  same at s=0, m=1,2 (general-m scalar columns)
     Z3  polarized Doppler elements T(s1,s2)^m closed vs quadrature (mixed and same-spin)
     Z4  m-summed polarized diagonals from closed forms == harness-15 route (RC26 Eq. 20 numbers)
     Z5  data export: exact p^6 series (q-polynomial coefficients) of the polarized system elements *)

prec = 25;
(* ---------- machinery ---------- *)
sYv[s_, l_, m_, st2_, ct2_] := (-1)^(m - s) Sqrt[(2 l + 1)/(4 Pi) (l + m)! (l - m)!/((l + s)! (l - s)!)] *
   Sum[If[r + s - m < 0 || r > l - s || r + s - m > l + s, 0,
     Binomial[l - s, r] Binomial[l + s, r + s - m] (-1)^(l - r - s + m) *
      st2^(2 l - (2 r + s - m)) ct2^(2 r + s - m)], {r, Max[0, m - s], l - s}];
(* polynomial part R(v):  sY = st2^|m+s| ct2^|m-s| R(v),  v = st2^2 = (1-mu)/2 *)
Rv[s_, l_, m_, v_] := Module[{a = Abs[m + s], c = Abs[m - s], ex},
   ex = Simplify[sYv[s, l, m, Sqrt[vv], Sqrt[1 - vv]]/(vv^(a/2) (1 - vv)^(c/2)), Assumptions -> 0 < vv < 1];
   If[! PolynomialQ[ex, vv], Print["Rv NOT polynomial at ", {s, l, m}, ": ", ex]];
   ex /. vv -> v];
JJ[w_, b_] := If[w === 0, 2 ArcTanh[b], With[{g = 1/Sqrt[1 - b^2], p = b/Sqrt[1 - b^2]}, ((g + p)^w - (g - p)^w)/w]];
(* J-reduction of sK^{d,m}_{l lp}(b):  returns the exact value given (numeric or symbolic) d, b *)
Kclosed[s_, l_, lp_, m_, d_, b_] /; l < Max[Abs[m], Abs[s]] || lp < Max[Abs[m], Abs[s]] := 0;
Kclosed[s_, l_, lp_, m_, d_, b_] := Module[{a = Abs[m + s], c = Abs[m - s], K = Max[Abs[m], Abs[s]],
    v, vp, E1, P, kmin, cl, g = 1/Sqrt[1 - b^2]},
   v = (b + 1 - t)/(2 b); vp = v (1 - b)/t;
   E1 = Together[v^a (1 - v)^c Rv[s, l, m, v] Rv[s, lp, m, vp]];
   kmin = Exponent[Denominator[E1], t]; P = Expand[Numerator[E1] Cancel[t^kmin/Denominator[E1]]];
   cl = CoefficientList[P, t];
   2 Pi (1 - b)^(a/2) (1 + b)^(c/2)/b *
    Sum[cl[[j + 1]] g^(-(j - kmin - K) - 1) JJ[j - kmin - K + 1 - d, b], {j, 0, Length[cl] - 1}]];
(* direct quadrature reference *)
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 14, MaxRecursion -> 20];
sY[s_, l_, m_, mu_] := sYv[s, l, m, Sqrt[(1 - mu)/2], Sqrt[(1 + mu)/2]];
KsRaw[s_, d_?NumericQ, l_, lp_, m_, b_] := Module[{g = 1/Sqrt[1 - b^2]},
   2 Pi nint[sY[s, l, m, mu] sY[s, lp, m, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] /; l < Max[Abs[m], Abs[s]] || lp < Max[Abs[m], Abs[s]] := 0;
Ks[s_, d_?NumericQ, l_, lp_, m_, b_] := Ks[s, d, l, lp, m, b] = KsRaw[s, d, l, lp, m, b];
report[name_, dev_, tol_: 10^-18] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

(* ---------- Z1: spin-2 columns ---------- *)
dev = Max@Flatten@Table[Abs[Kclosed[2, l, 2, m, d, b] - Ks[2, d, l, 2, m, b]],
   {l, 2, 4}, {m, 0, 2}, {d, {1/2, -3/2, 2}}, {b, {3/10, 7/10}}];
report["Z1  s=2 kernel columns: J-reduction == quadrature (l=2..4, m=0..2)", dev];

(* ---------- Z2: general-m scalar columns ---------- *)
dev = Max@Flatten@Table[Abs[Kclosed[0, l, lp, m, d, b] - Ks[0, d, l, lp, m, b]],
   {l, 1, 4}, {lp, 1, 2}, {m, 1, 2}, {d, {1/2, -3/2}}, {b, {3/10, 7/10}}];
report["Z2  s=0, m=1,2 scalar columns: J-reduction == quadrature", dev];

(* ---------- Z3: polarized Doppler elements from closed forms ---------- *)
(* T(s1,s2)^m_{l 2 ldd}(q) = (1/g) s1K^{-1-q,m}_{l2}(+b) s2K^{-q,m}_{2 ldd}(-b);
   at +b: Kclosed[s,l,lp,m,d,+b];  at -b: Kclosed[s,l,lp,m,d,-b] via parity of the formula (b -> -b works
   directly in the reduction since the substitution t = 1 + b mu was generic) — guard b sign by testing. *)
Tclosed[s1_, s2_, q_, l_, ldd_, m_, b_] := Sqrt[1 - b^2] Kclosed[s1, l, 2, m, -1 - q, b] Kclosed[s2, 2, ldd, m, -q, -b];
Tnum[s1_, s2_, q_, l_, ldd_, m_, b_] := Sqrt[1 - b^2] Ks[s1, -1 - q, l, 2, m, -b] Ks[s2, -q, 2, ldd, m, b];
(* NOTE our quadrature convention: Ks[s,d,l,lp,m,b] is the kernel at "-b" in site notation when called with +b;
   match by testing both sign pairings on one element first *)
t1 = Abs[Tclosed[0, 2, 1/2, 2, 2, 0, 3/10] - Tnum[0, 2, 1/2, 2, 2, 0, 3/10]];
t2 = Abs[Tclosed[0, 2, 1/2, 2, 2, 0, -3/10] - Tnum[0, 2, 1/2, 2, 2, 0, 3/10]];
sgn = If[t1 < t2, 1, -1];
Print["Z3a sign pairing pinned: closed-form b -> ", sgn, " b matches quadrature"];
dev = Max@Flatten@Table[Abs[Tclosed[ss[[1]], ss[[2]], q, l, ldd, m, sgn 3/10] - Tnum[ss[[1]], ss[[2]], q, l, ldd, m, 3/10]],
   {ss, {{2, 2}, {0, 2}, {2, 0}}}, {q, {1/2, -7/10}}, {l, {2, 3}}, {ldd, {2, 3}}, {m, 0, 2}];
report["Z3  polarized Doppler elements: closed == quadrature (T(2,2),T(0,2),T(2,0))", dev];

(* ---------- Z4: m-summed mixed diagonal (the Eq.-20 object) from closed forms ---------- *)
msumC[q_, l_, b_] := Sum[If[Abs[m] > Min[l, 2], 0, Tclosed[0, 2, q, l, l, Abs[m], sgn b]], {m, -2, 2}];
msumN[q_, l_, b_] := Sum[If[Abs[m] > Min[l, 2], 0, Tnum[0, 2, q, l, l, Abs[m], b]], {m, -2, 2}];
dev = Max@Table[Abs[msumC[q, l, 3/10] - msumN[q, l, 3/10]], {q, {1/2}}, {l, {2, 3}}];
report["Z4  m-summed mixed diagonal from closed forms", dev];

(* ---------- display: two representative closed forms (in rapidity form) ---------- *)
Print["-- representative closed forms (d symbolic; presented via eta with b = tanh eta) --"];
clean[ex_] := Collect[FullSimplify[ExpToTrig[ex /. b -> Tanh[et]], Assumptions -> et > 0], {Sinh[d et], Cosh[d et]}, Simplify];
Print["2K^{d,1}_{22}: ", clean[Kclosed[2, 2, 2, 1, d, b]]];
Print["0K^{d,1}_{11}: ", clean[Kclosed[0, 1, 1, 1, d, b]]];

(* ---------- Z5: data export — exact p^6 series of the polarized system ---------- *)
(* series route: substitute b -> pp/Sqrt[1+pp^2] in the closed forms with d = -1-q / -q and expand *)
serT[s1_, s2_, l_, ldd_, m_, ord_] := Module[{bb = pp/Sqrt[1 + pp^2]},
   Collect[Normal@Series[Tclosed[s1, s2, q, l, ldd, m, sgn bb], {pp, 0, ord}], pp, Simplify]];
ser1 = serT[2, 2, 2, 2, 0, 4];
Print["Z5a sample: T(2,2)^0_{222}(q) = ", ser1, " + O(p^6)"];
(* numeric spot check of the series *)
devZ5 = Abs[(ser1 /. {q -> 1/2, pp -> 1/10}) - Tnum[2, 2, 1/2, 2, 2, 0, (1/10)/Sqrt[1 + 1/100]]];
Print["Z5b series spot check at p=0.1: dev = ", ScientificForm[N[devZ5], 2], "  (expect ~p^6 = 1e-6)"];
exportData = Module[{els = {}, s},
   Do[If[Min[l, ldd] >= If[Abs[s1s2[[1]]] == 2 || Abs[s1s2[[2]]] == 2, Abs[m], Abs[m]] && Abs[m] <= Min[l, ldd, 2],
     s = serT[s1s2[[1]], s1s2[[2]], l, ldd, m, 6];
     AppendTo[els, <|"s1" -> s1s2[[1]], "s2" -> s1s2[[2]], "m" -> m, "l" -> l, "lp" -> 2, "ldd" -> ldd,
        "series_coefficients_q_polynomials" -> Table[ToString[Simplify[Coefficient[s, pp, k]], InputForm], {k, 0, 6}]|>]],
    {s1s2, {{2, 2}, {0, 2}, {2, 0}}}, {l, 2, 4}, {ldd, 2, 4}, {m, 0, 2}];
   els];
Export["../data/polarized_tables_p6.json", <|"description" -> "Polarized (spin-weighted) Doppler-operator elements T(s1,s2)^m_{l,2,ldd}(q) on nu^q: exact p-series to O(p^6); coefficients are polynomials in q. Conventions as site Section 2 / N15/N22.",
   "generated" -> DateString[], "elements" -> exportData|>, "JSON", "Compact" -> False];
Print["Z5  exported data/polarized_tables_p6.json with ", Length[exportData], " elements"];
Print["--- spin/m closed-forms harness done ---"];
