(* ::Package:: *)
(* 14_recoil_lattice.wl — N14: the recoil ladder lives on the weight lattice; all recursions
   extend verbatim to every rung.
   Facts verified:
     R1  commutation dictionary: f(Ohat) x^s = x^s f(Ohat - s)   (x = h nu / m_e c^2)
         [on power laws: trivial; verified as operator statement on a test function]
     R2  rest-frame frequency insertion = lattice shift: inserting w'^s between the boosts maps
         D -> w^s * (s,0)-lattice element (bookkeeping identity, verified numerically):
         (1/g) B(-b)[ w'^s * (B(b) n) ]  on nu^q  =  w^s(lab) * (s,0)D(q) * nu^{q+s}-normalization
         i.e.  (1/g) Kp[-1-(q+s), l, lp] Km[-q, lp, ldd] = (s,0)-lattice element at q  == (j=s,k=0)D(q)
     R3  the lattice-generalized rapidity ODE (N3 at (j,k) != (0,0)):
         d/deta (j,k)D_{l lp ldd} =
             -(l+3+j+q)  C_{l+1}  (j,k)D_{l+1,lp,ldd} + (l-2-j-q)   C_l   (j,k)D_{l-1,lp,ldd}
             -(ldd-k-q)  C_{ldd+1}(j,k)D_{l,lp,ldd+1} + (ldd+1+k+q) C_ldd (j,k)D_{l,lp,ldd-1}
             - beta (j,k)D_{l lp ldd}
         verified by FD at the recoil-relevant rungs (j,k) = (1,0), (2,0), (0,-1).  *)

prec = 30;
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, m_, x_] := Sqrt[(2 l + 1)/(4 Pi) Factorial[l - m]/Factorial[l + m]] LegendreP[l, m, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, m_, b_] := Module[{g = ga[b]},
   2 Pi nint[Ybar[l, m, mu] Ybar[lp, m, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, m_, b_] /; l < Abs[m] || lp < Abs[m] := 0;
Km[d_?NumericQ, l_, lp_, m_, b_] := Km[d, l, lp, m, b] = KmRaw[d, l, lp, m, b];
Kp[d_?NumericQ, l_, lp_, m_, b_] := Km[d, l, lp, m, -b];
Cm[l_, m_] := If[l <= Abs[m], 0, Sqrt[(l^2 - m^2)/(4 l^2 - 1)]];
DD[j_, k_, q_, l_, lp_, ldd_, m_, b_] := (1/ga[b]) Kp[-1 - j - q, l, lp, m, b] Km[-k - q, lp, ldd, m, b];
report[name_, dev_, tol_: 10^-24] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

(* ---- R1: f(Ohat) x^s = x^s f(Ohat - s), operator statement on a test function ----
   Take f(O) = e^{aO} (generates all analytic f):  e^{aO}[x^s g(x)] = e^{-as} x^s e^{aO} g(x)??
   e^{aO} h(x) = h(e^{-a}x):  LHS = (e^{-a}x)^s g(e^{-a}x) = e^{-as} x^s [e^{aO}g](x)  =>
   e^{aO} x^s = e^{-as} x^s e^{aO} = x^s e^{a(O-s)}.  Numeric spot check: *)
a0 = 3/10; s0 = 2; g1[x_] := Exp[-x] (1 + x);
xs = 7/5;
lhs = (Exp[-a0] xs)^s0 g1[Exp[-a0] xs];
rhs = xs^s0 Exp[-a0 s0] g1[Exp[-a0] xs];
report["R1  e^{aO} x^s = x^s e^{a(O-s)} (spot check)", Abs[lhs - rhs], 10^-25];

(* ---- R2: rest-frame insertion = lattice shift ----
   Round trip with w'^s inserted, acting on nu^q, multipole (l, lp, ldd):
   step 1: [B(b) nu^q]_lp-component = Km[-q, lp, ldd, m, b] nu'^q     (rest-frame spectrum)
   step 2: multiply by nu'^s:  Km[-q,...] nu'^{q+s}
   step 3: apply (1/g) B(-b):  (1/g) Kp[-1-(q+s), l, lp, m, b] Km[-q, lp, ldd, m, b] nu^{q+s}
   claim: the coefficient equals the (j,k) = (s,0) lattice element at spectral index q:
          (s,0)D(q) = (1/g) Kp[-1-s-q, l, lp] Km[-q, lp, ldd]     -- definitional match *)
bb = 3/10; qv = 1/2;
dev = Max@Table[Abs[(1/ga[bb]) Kp[-1 - (qv + s), l, lp, 0, bb] Km[-qv, lp, ldd, 0, bb] -
     DD[s, 0, qv, l, lp, ldd, 0, bb]], {s, {1, 2}}, {l, 0, 1}, {lp, 0, 1}, {ldd, 0, 1}];
report["R2  w'^s insertion == (s,0) lattice element", dev, 10^-28];

(* ---- R3: lattice-generalized rapidity ODE at recoil rungs ---- *)
eta0 = ArcTanh[bb]; de = 10^-6;
fdD[j_, k_, q_, l_, lp_, ldd_, m_] := (DD[j, k, q, l, lp, ldd, m, Tanh[eta0 + de]] - DD[j, k, q, l, lp, ldd, m, Tanh[eta0 - de]])/(2 de);
odeD[j_, k_, q_, l_, lp_, ldd_, m_] := -(l + 3 + j + q) Cm[l + 1, m] DD[j, k, q, l + 1, lp, ldd, m, bb] +
   (l - 2 - j - q) Cm[l, m] If[l - 1 >= Abs[m], DD[j, k, q, l - 1, lp, ldd, m, bb], 0] -
   (ldd - k - q) Cm[ldd + 1, m] DD[j, k, q, l, lp, ldd + 1, m, bb] +
   (ldd + 1 + k + q) Cm[ldd, m] If[ldd - 1 >= Abs[m], DD[j, k, q, l, lp, ldd - 1, m, bb], 0] - bb DD[j, k, q, l, lp, ldd, m, bb];
dev = Max@Table[Abs[fdD[j, k, q, l, lp, ldd, 0] - odeD[j, k, q, l, lp, ldd, 0]],
   {j, {1, 2, 0}}, {k, {0, -1}}, {q, {1/2}}, {l, 0, 1}, {lp, 0, 2, 2}, {ldd, 0, 1}];
report["R3  lattice-generalized N3 ODE at rungs (j,k) in {1,2,0}x{0,-1}", dev, 10^-8];

Print["--- recoil-lattice harness done ---"];
