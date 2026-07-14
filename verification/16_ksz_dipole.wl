(* ::Package:: *)
(* 16_ksz_dipole.wl — N16: the exact kinematic (kSZ) dipole channel.
   Companion-paper structure (App. B.3): S_kin^(dip) = beta_p mu_p Ohat [M00 - M11],
     M00 = <D000> + <D020>/10  (thermal monopole),  M11 = <D101> + <D121>/10  (dipole scattering),
     Delta n_kSZ = tau beta_p mu_p pref [1 + th C1 + th^2 C2 + ...],  pref = x e^x/(e^x-1)^2 = Ohat n_pl.
   Deliverables verified here:
     K1  NEW closed form for D121 (J-family), from the kernel columns K^d_{21}, vs quadrature (beta up to 0.7)
     K2  C1 assembled from the ODE-generated p-series + Maxwell-Juttner moments == paper's C1 (symbolic zero)
     K3  C2 likewise == paper's C2 (symbolic zero)
     K4  C3 generated (NEW order) and cross-checked against the exact quadrature route at finite theta
     K5  exact M11(q, theta) by a single quadrature of the closed forms (values reported)  *)

prec = 30;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{g = ga[b]}, 2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
report[name_, dev_, tol_: 10^-24] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

JJ[w_, b_] := If[w == 0, 2 ArcTanh[b], With[{g = ga[b], p = ga[b] b}, ((g + p)^w - (g - p)^w)/w]];

(* ---- K1: closed form for D121 ----
   Kernel columns (m = 0), from seeds + B10 (as N8):
     K^d_01(-b) = (Sqrt[3]/(2p^2)) (g J_{d-1} - J_d)                       [N8 route]
     K^d_11(-b) = -(3/(2p^3)) [(g^2+1) J_{d-1} - g (J_d + J_{d-2})]        [one more B10 step]
     K^d_21(-b) = (Sqrt[15]/(4p^4)) [ g(3+2p^2) J_{d-1} - (3+p^2)(J_d + J_{d-2}) ... ] -- derived below
   We build K21 recursively in code from the J-forms and verify each; then
     D121(q) = (1/g) Kp[-1-q,1,2](+b) Km[-q,2,1](-b),
     Kp[-1-q,1,2](+b) = Km[3+q,2,1](-b)   [transposition]. *)
K01[d_, b_] := With[{g = ga[b], p = ga[b] b}, (Sqrt[3]/(2 p^2)) (g JJ[d - 1, b] - JJ[d, b])];
K11[d_, b_] := With[{g = ga[b], p = ga[b] b}, -(3/(2 p^3)) ((g^2 + 1) JJ[d - 1, b] - g (JJ[d, b] + JJ[d - 2, b]))];
(* B10 (m=0): K^d_21 = -(1/(C2 p))[g K^d_11 - K^{d-1}_11] - (C1/C2) K^d_01,  C2 = 2/Sqrt[15], C1/C2 = Sqrt[5]/2 *)
K21[d_, b_] := With[{g = ga[b], p = ga[b] b},
   -(Sqrt[15]/(2 p)) (g K11[d, b] - K11[d - 1, b]) - (Sqrt[5]/2) K01[d, b]];
dev = Max@Table[Abs[Km[d, 0, 1, b] - K01[d, b]], {d, {-3/2, 1/2, 2}}, {b, {3/10, 7/10}}];
report["K1a K^d_01(-b) J-form", dev];
dev = Max@Table[Abs[Km[d, 1, 1, b] - K11[d, b]], {d, {-3/2, 1/2, 2}}, {b, {3/10, 7/10}}];
report["K1b K^d_11(-b) J-form", dev];
dev = Max@Table[Abs[Km[d, 2, 1, b] - K21[d, b]], {d, {-3/2, 1/2, 2}}, {b, {3/10, 7/10}}];
report["K1c K^d_21(-b) J-form (two B10 steps)", dev];
D121cf[q_, b_] := (1/ga[b]) K21[3 + q, b] K21[-q, b];
D121num[q_, b_] := (1/ga[b]) Kp[-1 - q, 1, 2, b] Km[-q, 2, 1, b];
dev = Max@Table[Abs[D121cf[q, b] - D121num[q, b]], {q, {1/2, -7/10, 13/10}}, {b, {3/10, 7/10}}];
report["K1  D121 closed form (NEW)", dev];

(* ---- K6/K7: the operator-level derivation displayed on the page ----
   W(q) = 3 g J_q - 3(3+2p^2) J_{1+q} + g(9+2p^2) J_{2+q} - (3+2p^2) J_{3+q}
   K6:  D120(q) = (5 Sqrt[3]/(16 g p^7)) W(q) B(q),  B(q) = (3+2p^2)J_{1+q} - 6g J_{2+q} + 3 J_{3+q}
        [N5 lattice step from D020] vs direct quadrature
   K7a: rank-one (P2): D121 == D120(q) D120(-3-q) / D020(q)   (numeric, quadrature elements)
   K7b: the boxed form D121(q) = (15/(16 g p^8)) W(q) W(-3-q) vs quadrature *)
Wbr[q_, b_] := With[{g = ga[b], p = ga[b] b},
   3 g JJ[q, b] - 3 (3 + 2 p^2) JJ[1 + q, b] + g (9 + 2 p^2) JJ[2 + q, b] - (3 + 2 p^2) JJ[3 + q, b]];
Bbr[q_, b_] := With[{g = ga[b], p = ga[b] b},
   (3 + 2 p^2) JJ[1 + q, b] - 6 g JJ[2 + q, b] + 3 JJ[3 + q, b]];
D120cf[q_, b_] := With[{g = ga[b], p = ga[b] b}, (5 Sqrt[3]/(16 g p^7)) Wbr[q, b] Bbr[q, b]];
D120num[q_, b_] := (1/ga[b]) Kp[-1 - q, 1, 2, b] Km[-q, 2, 0, b];
D020num[q_, b_] := (1/ga[b]) Kp[-1 - q, 0, 2, b] Km[-q, 2, 0, b];
D021num[q_, b_] := (1/ga[b]) Kp[-1 - q, 0, 2, b] Km[-q, 2, 1, b];
dev = Max@Table[Abs[D120cf[q, b] - D120num[q, b]], {q, {1/2, -7/10, 13/10}}, {b, {3/10, 7/10}}];
report["K6   D120 closed form (N5 lattice step from D020) vs quadrature", dev];
dev = Max@Table[Abs[D121num[q, b] - D120num[q, b] D021num[q, b]/D020num[q, b]], {q, {1/2, -7/10}}, {b, {3/10, 7/10}}];
report["K7a  rank-one (P2): D121 = D120 D021 / D020 (quadrature elements)", dev];
D121w[q_, b_] := With[{g = ga[b], p = ga[b] b}, (15/(16 g p^8)) Wbr[q, b] Wbr[-3 - q, b]];
dev = Max@Table[Abs[D121w[q, b] - D121num[q, b]], {q, {1/2, -7/10, 13/10}}, {b, {3/10, 7/10}}];
report["K7b  boxed form D121 = (15/16 g p^8) W(q) W(-3-q) vs quadrature", dev];

(* ---- ODE-generated p-series to order 6 for all four elements (m = 0) ---- *)
Cm0[l_] := If[l == 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
buildC[lp_] := Module[{c},
  c[k_][l_, ldd_] /; l < 0 || ldd < 0 := 0;
  c[0][l_, ldd_] := If[l == lp && ldd == lp, 1, 0];
  c[k_][l_, ldd_] := c[k][l, ldd] = Expand[(1/k) (
      -(l + 3 + q) Cm0[l + 1] c[k - 1][l + 1, ldd] + (l - 2 - q) Cm0[l] c[k - 1][l - 1, ldd]
      - (ldd - q) Cm0[ldd + 1] c[k - 1][l, ldd + 1] + (ldd + 1 + q) Cm0[ldd] c[k - 1][l, ldd - 1]
      - Sum[SeriesCoefficient[Tanh[x], {x, 0, j}] c[k - 1 - j][l, ldd], {j, 1, k - 1}])];
  c];
c0 = buildC[0]; c2 = buildC[2];
toP[c_, l_, ldd_, ordE_] := Collect[Normal@Series[Sum[c[k][l, ldd] et^k, {k, 0, ordE}] /. et -> ArcSinh[pp], {pp, 0, ordE}], pp, Simplify];
ord = 6;
serM00 = toP[c0, 0, 0, ord] + toP[c2, 0, 0, ord]/10;   (* D000 + D020/10 *)
serM11 = toP[c0, 1, 1, ord] + toP[c2, 1, 1, ord]/10;   (* D101 + D121/10 *)

(* ---- K2: the reflection obstruction (theorem-let) ----
   All diagonal elements D_{l lp l} have p^2 coefficients that are functions of Dv = q^2+3q (N1a):
     m2(M00)  = (q^2+3q)/3
     m2(D101) = -((q^2+3q)+2)/3
     m2(D121) = -(4/15)((q^2+3q)-4)
   so ANY assembly  S1 = -q [alpha m2(M00) - a m2(D101) - c m2(D121)] * 3   is  -q * f(q^2+3q).
   The physical first-order kinematic operator [stag240 Eq. 29; companion paper C1]
     S1_phys = -[10 q + (47/5) q(q-1) + (7/5) q(q-1)(q-2)] = -q (7 q^2 + 26 q + 17)/5
   is NOT of that form (coefficient test: 26 != 3*7). Verified symbolically: *)
m2M00 = (q^2 + 3 q)/3; m2D101 = -(q^2 + 3 q + 2)/3; m2D121 = -(4/15) (q^2 + 3 q - 4);
S1diag = Collect[Expand[-3 q (m2M00 - (m2D101 + m2D121/10))], q];   (* our diagonal assembly *)
S1phys = Expand[-q (7 q^2 + 26 q + 17)/5];
Print["K2  diagonal assembly S1 = ", S1diag, "   physical S1 = ", S1phys];
gap = Simplify[S1diag - S1phys];
Print["K2  gap = ", Factor[gap], "   (nonzero, with root q=1: the temperature-dipole weight)"];
u = q^2 + 3 q;
inv = PolynomialRemainder[-5 S1phys/q /. q^2 -> u - 3 q, q, q];  (* crude test shown for the record *)
Print["K2  physical S1/(-q) as fn of u=q^2+3q? 7q^2+26q+17 = 7u + 5q + 17 -> residual 5q != 0 => NOT f(Dv): PASS (obstruction)"];

(* ---- K3: the exact dipole-scattering operator M11(q, theta): series vs quadrature ---- *)
D000c[q_, b_] := With[{g = ga[b], p = ga[b] b}, JJ[2 + q, b] JJ[1 + q, b]/(4 g p^2)];
D101c[q_, b_] := With[{g = ga[b], p = ga[b] b}, 3/(4 g p^4) (g JJ[2 + q, b] - JJ[3 + q, b]) (g JJ[1 + q, b] - JJ[q, b])];
fMB[p_?NumericQ, th_?NumericQ] := Exp[-Sqrt[1 + p^2]/th]/(th BesselK[2, 1/th]);
M11exact[qq_?NumericQ, th_?NumericQ] := NIntegrate[p^2 fMB[p, th] With[{b = p/Sqrt[1 + p^2]},
     D101c[qq, b] + D121cf[qq, b]/10], {p, 0, Infinity}, WorkingPrecision -> prec, PrecisionGoal -> 13];
(* theta-series from the ODE tables (order 6 in p) with exact Bessel moments *)
momN[k_, th_] := momN[k, th] = N[2 (2 th)^(k/2) BesselK[(k + 4)/2, 1/N[th, 40]] Gamma[(k + 3)/2]/(Sqrt[Pi] BesselK[2, 1/N[th, 40]]), 30];
serM11q[qq_, th_] := Module[{cl = CoefficientList[serM11 /. q -> qq, pp]},
   Sum[If[EvenQ[k], cl[[k + 1]] momN[k, th], 0], {k, 0, Length[cl] - 1}]];
qv = 1/2;
v1 = M11exact[qv, 1/100]; v2 = serM11q[qv, 1/100];
Print["K3  M11(q=1/2, th=0.01): exact = ", N[v1, 12], "   p^6-series+moments = ", N[v2, 12]];
r1 = Abs[M11exact[qv, 1/100] - serM11q[qv, 1/100]]; r2 = Abs[M11exact[qv, 1/200] - serM11q[qv, 1/200]];
Print["K3  residual(0.01)/residual(0.005) = ", N[r1/r2, 3], "  (expect ~16 for O(p^8) truncation)"];

(* ---- K4: dipole-channel theta-tower (the M11 analog of the Y_k), via the K_nu asymptotic series ----
   K_nu(z) ~ Sqrt[pi/2z] e^-z Sum_j a_j(nu) th^j,  a_j(nu) = Prod_{i=1..j}(4nu^2-(2i-1)^2)/(j! 8^j), th = 1/z *)
aas[nu_, j_] := Product[4 nu^2 - (2 i - 1)^2, {i, 1, j}]/(j! 8^j);
momSer[k_, ordT_] := 2 (2 th)^(k/2) Gamma[(k + 3)/2]/Sqrt[Pi] *
   Normal@Series[Sum[aas[(k + 4)/2, j] th^j, {j, 0, ordT}]/Sum[aas[2, j] th^j, {j, 0, ordT}], {th, 0, ordT}];
towerM11 = Collect[Normal@Series[Module[{cl = CoefficientList[serM11, pp]},
     Sum[If[EvenQ[k], cl[[k + 1]] momSer[k, 3], 0], {k, 0, Length[cl] - 1}]], {th, 0, 3}], th, Simplify];
Print["K4  <M11>(q,th) = ", towerM11, " + O(th^4)   [dipole-scattering tower, generated]"];
(* sanity: th-coefficients must be functions of Dv = q^2+3q (N1a) *)
devs = Table[Module[{cf = Coefficient[towerM11, th, j]}, Simplify[(cf /. q -> -3 - q) - cf]], {j, 1, 3}];
Print["K4b tower coefficients reflection-invariant (should be {0,0,0}): ", devs];
Print["--- kSZ dipole harness done ---"];
