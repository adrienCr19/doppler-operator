(* ::Package:: *)
(* 02_doppler.wl — Verify NEW recurrence relations & symmetries for the Doppler operator
      D^m_{l lp ldd}(nu, beta) = (1/gamma) ^{-1}B^m_{l lp}(nu,-beta) ^{0}B^m_{lp ldd}(nu,beta)
   [Rosenberg & Chluba 2026 (stag331) Eq. 9; Chluba & Rosenberg 2026 (stag240) Eq. 3]
   On power-law test spectra nu^q (eigenfunctions: Ohat nu^q = -q nu^q) all operator weights
   become scalars, so operator identities reduce to kernel-product identities:
      D(q)_{l lp ldd} = (1/g) Kp[-1-q, l, lp] Km[-q, lp, ldd],
   Kp/Km = aberration kernel at +beta/-beta as in 01_conventions.wl.  s = 0 file. *)

prec = 30;
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, m_, x_] := Sqrt[(2 l + 1)/(4 Pi) Factorial[l - m]/Factorial[l + m]] LegendreP[l, m, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, m_, b_] := Module[{g = ga[b]},
   2 Pi nint[ Ybar[l, m, mu] Ybar[lp, m, (mu + b)/(1 + b mu)] / (g (1 + b mu))^d ]];
Km[d_?NumericQ, l_, lp_, m_, b_] /; l < Abs[m] || lp < Abs[m] := 0;
Km[d_?NumericQ, l_, lp_, m_, b_] := Km[d, l, lp, m, b] = KmRaw[d, l, lp, m, b];
Kp[d_?NumericQ, l_, lp_, m_, b_] := Km[d, l, lp, m, -b];
Cm[l_, m_] := If[l <= Abs[m], 0, Sqrt[(l^2 - m^2)/(4 l^2 - 1)]];
report[name_, dev_, tol_: 10^-10] := Print[name, ": ", If[dev < tol, "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

(* Doppler operator on power law nu^q *)
DD[q_, l_, lp_, ldd_, m_, b_] := (1/ga[b]) Kp[-1 - q, l, lp, m, b] Km[-q, lp, ldd, m, b];

bb = 3/10; qs = {-7/10, 1/2, 13/10};

(* ---- D1: corrected seed  ^dK^m_{mm}(-b) = g^{-(m+d)} 2F1((m+d)/2,(m+d+1)/2; 3/2+m; b^2)  [stag331 B4] ---- *)
dev = Max@Table[Abs[Km[d, m, m, m, bb] - ga[bb]^-(m + d) Hypergeometric2F1[(m + d)/2, (m + d + 1)/2, 3/2 + m, bb^2]],
   {d, {-1, 0, 1, 2, 7/3}}, {m, 0, 3}];
report["D1  seed K^m_mm 2F1 (stag331 B4, corrected reading)", dev];

(* ---- D2: transposition in this convention:  ^dK_{l lp}(-b) = ^{2-d}K_{lp l}(+b)  [equiv. stag698 13a] ---- *)
dev = Max@Table[Abs[Km[d, l, lp, m, bb] - Kp[2 - d, lp, l, m, bb]], {d, {0, 1, 2, -1, 1/2}}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}];
report["D2  kernel transposition K^d_{llp}(-b)=K^{2-d}_{lpl}(b)", dev];

(* ---- D3: column-index weight relation (corrected HLC26 A.7b):
   Kp[d+1, lp, l] = g Kp[d, lp, l] + g b (C_{l+1} Kp[d, lp, l+1] + C_l Kp[d, lp, l-1]) ---- *)
dev = Max@Table[Abs[Kp[d + 1, lp, l, m, bb] - ga[bb] (Kp[d, lp, l, m, bb] + bb (Cm[l + 1, m] Kp[d, lp, l + 1, m, bb] + Cm[l, m] If[l - 1 >= Abs[m], Kp[d, lp, l - 1, m, bb], 0]))],
   {d, {0, 1, -1, 3/2}}, {m, 0, 1}, {lp, m, 3}, {l, m, 3}];
report["D3  weight raising, column-shift form", dev];

(* ============ NEW RESULTS: Doppler operator ============ *)

(* ---- N1: REFLECTION SYMMETRY   D_{l lp ldd}(q) = D_{ldd lp l}(-3-q)  i.e. Ohat -> 3-Ohat ---- *)
dev = Max@Table[Abs[DD[q, l, lp, ldd, m, bb] - DD[-3 - q, ldd, lp, l, m, bb]],
   {q, qs}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}, {ldd, m, 3}];
report["N1  reflection symmetry D(Ohat) = D^T(3-Ohat)", dev];

(* ---- N2: PARITY   D_{l lp ldd}(-beta) = (-1)^(l+ldd) D_{l lp ldd}(beta) ---- *)
dev = Max@Table[Abs[DD[q, l, lp, ldd, m, -bb] - (-1)^(l + ldd) DD[q, l, lp, ldd, m, bb]],
   {q, {1/2}}, {m, 0, 1}, {l, m, 3}, {lp, m, 3}, {ldd, m, 3}];
report["N2  parity D(-beta) = (-1)^(l+l'') D(beta)", dev];

(* ---- N3: RAPIDITY ODE (closed at fixed intermediate lp):
   d/deta D_{l lp ldd} = -(l+3+q) C_{l+1} D_{l+1,lp,ldd} + (l-2-q) C_l D_{l-1,lp,ldd}
                         -(ldd-q) C_{ldd+1} D_{l,lp,ldd+1} + (ldd+1+q) C_{ldd} D_{l,lp,ldd-1}
                         - beta D_{l lp ldd}
   (operator form: q -> -Ohat).  Verified against central finite differences in eta. ---- *)
eta0 = ArcTanh[bb]; de = 10^-6;
fdD[q_, l_, lp_, ldd_, m_] := (DD[q, l, lp, ldd, m, Tanh[eta0 + de]] - DD[q, l, lp, ldd, m, Tanh[eta0 - de]])/(2 de);
odeD[q_, l_, lp_, ldd_, m_] := -(l + 3 + q) Cm[l + 1, m] DD[q, l + 1, lp, ldd, m, bb] +
   (l - 2 - q) Cm[l, m] If[l - 1 >= Abs[m], DD[q, l - 1, lp, ldd, m, bb], 0] -
   (ldd - q) Cm[ldd + 1, m] DD[q, l, lp, ldd + 1, m, bb] +
   (ldd + 1 + q) Cm[ldd, m] If[ldd - 1 >= Abs[m], DD[q, l, lp, ldd - 1, m, bb], 0] - bb DD[q, l, lp, ldd, m, bb];
dev = Max@Table[Abs[fdD[q, l, lp, ldd, m] - odeD[q, l, lp, ldd, m]],
   {q, {1/2, -7/10}}, {m, 0, 1}, {l, m, 2}, {lp, m, 2}, {ldd, m, 2}];
report["N3  rapidity ODE for Doppler operator (FD check)", dev, 10^-8];

(* ---- N4: SUM RULE (cf. stag698 Sec 4.2 / stag331 Eq. 10-11):
   Sum_lp D_{l lp ldd} = delta_{l,ldd} - beta (C^m_{ldd+1} delta_{l,ldd+1} + C^m_{ldd} delta_{l,ldd-1}) ---- *)
Lmax = 14;
dev = Max@Table[Abs[Sum[DD[q, l, lp, ldd, m, bb], {lp, Abs[m], Lmax}] -
     (KroneckerDelta[l, ldd] - bb (Cm[ldd + 1, m] KroneckerDelta[l, ldd + 1] + Cm[ldd, m] KroneckerDelta[l, ldd - 1]))],
   {q, {1/2}}, {m, 0, 1}, {l, m, 2}, {ldd, m, 2}];
Print["N4  sum rule over intermediate lp, Lmax=", Lmax, ": max dev = ", ScientificForm[N[dev], 3], " (truncation-limited)"];

(* ---- N5: WEIGHT-LATTICE l-RAISING (outgoing index; companion operator has first weight -2):
   Define companion  E(q)_{l lp ldd} = (1/g) Kp[-2-q, l, lp] Km[-q, lp, ldd]   (= (1/g) ^{-2}B(-b) ^0B(b) el.)
   Then from kernel recursion (stag331 B3 + parity):
   D_{l lp ldd} = (1/(C_l p)) ( g D_{l-1,lp,ldd} - E_{l-1,lp,ldd} ) - (C_{l-1}/C_l) D_{l-2,lp,ldd} ---- *)
EE[q_, l_, lp_, ldd_, m_, b_] := (1/ga[b]) Kp[-2 - q, l, lp, m, b] Km[-q, lp, ldd, m, b];
dev = Max@Table[With[{p = ga[bb] bb},
    Abs[DD[q, l, lp, ldd, m, bb] - ((1/(Cm[l, m] p)) (ga[bb] DD[q, l - 1, lp, ldd, m, bb] - EE[q, l - 1, lp, ldd, m, bb]) -
        (Cm[l - 1, m]/Cm[l, m]) If[l - 2 >= Abs[m], DD[q, l - 2, lp, ldd, m, bb], 0])]],
   {q, {1/2, -7/10}}, {m, 0, 1}, {l, Max[Abs[m], 1] + 1, 3}, {lp, m, 2}, {ldd, m, 2}];
report["N5  outgoing-index raising w/ weight companion", dev];

(* ---- N6: incoming-index (ldd) raising via reflection of N5:
   D_{l lp ldd}(q) = (1/(C_ldd p)) ( g D_{l,lp,ldd-1}(q) - F_{l,lp,ldd-1}(q) ) - (C_{ldd-1}/C_ldd) D_{l,lp,ldd-2}(q)
   with F(q)_{l lp ldd} = (1/g) Kp[-1-q, l, lp] Km[-q+1, lp, ldd]  (second weight lowered by 1: reflection image of E) ---- *)
FF[q_, l_, lp_, ldd_, m_, b_] := (1/ga[b]) Kp[-1 - q, l, lp, m, b] Km[1 - q, lp, ldd, m, b];
dev = Max@Table[With[{p = ga[bb] bb},
    Abs[DD[q, l, lp, ldd, m, bb] - ((1/(Cm[ldd, m] p)) (ga[bb] DD[q, l, lp, ldd - 1, m, bb] - FF[q, l, lp, ldd - 1, m, bb]) -
        (Cm[ldd - 1, m]/Cm[ldd, m]) If[ldd - 2 >= Abs[m], DD[q, l, lp, ldd - 2, m, bb], 0])]],
   {q, {1/2, -7/10}}, {m, 0, 1}, {ldd, Max[Abs[m], 1] + 1, 3}, {lp, m, 2}, {l, m, 2}];
report["N6  incoming-index raising w/ weight companion", dev];

(* ---- N7: companion reduction: E is itself a Doppler-type operator; check E(q) = gamma-relation via D3:
   Kp[-2-q, l, lp] = g( Kp[-1-q,l,lp] - b(C_{lp+1}Kp[-1-q,l,lp+1] + C_lp Kp[-1-q,l,lp-1]) )  [R1 row->col?]
   Using D3 with d+1 = -1-q at column index lp:
   Kp[-1-q, l, lp] = g Kp[-2-q, l, lp] + g b (C_{lp+1} Kp[-2-q, l, lp+1] + C_lp Kp[-2-q, l, lp-1])
   => E couples back to D with intermediate-index shifts:
   D_{l lp ldd} = g E_{l lp ldd} + g b ( C_{lp+1} (1/g)Kp[-2-q,l,lp+1]Km[-q,lp,ldd] + ... )   [mismatched: NOT closed]
   Instead verify the *correct* closed statement: E_{l lp ldd}(q) = D_{l lp ldd}(q+1) * scaling?
   ^{-2}B(nu,-b) nu^q = Kp[-2-q] nu^q and ^{-1}B(nu,-b) nu^{q+1} = Kp[-2-q] nu^{q+1}:
   so E(q)_{l lp ldd} = (1/g) Kp[-1-(q+1), l, lp] Km[-q, lp, ldd]
                      = D(q+1)_{l lp ldd} * Km[-q,lp,ldd]/Km[-q-1,lp,ldd]   (elementwise, NOT operator id)
   -> No closed elementwise identity expected; we verify instead the OPERATOR statement:
      E(q) = (1/g) Kp[-2-q,l,lp] Km[-q,lp,ldd]  equals  the (j,k)=(1,0) lattice point, def-check only. ---- *)
Print["N7  (companion is lattice point (j,k)=(1,0); no elementwise reduction — by design)"];

(* ---- N8: symbolic Taylor tables from the ODE, checked against direct numerics.
   Generate D(q)_{l lp ldd}(eta) as series in eta to O(eta^4) using N3 ODE with beta = tanh(eta),
   seed D|_{eta=0} = delta_{l lp} delta_{lp ldd}.  Then compare at beta=0.05 against direct integration. ---- *)
LmaxS = 8;
buildSeries2[lp_, q_, ordE_] := Module[{c},
  c[k_][l_, ldd_] /; l < 0 || ldd < 0 := 0;
  c[0][l_, ldd_] := If[l == lp && ldd == lp, 1, 0];
  c[k_][l_, ldd_] := c[k][l, ldd] = Simplify[(1/k) (
      -(l + 3 + q) Cm[l + 1, 0] c[k - 1][l + 1, ldd] + (l - 2 - q) Cm[l, 0] c[k - 1][l - 1, ldd]
      - (ldd - q) Cm[ldd + 1, 0] c[k - 1][l, ldd + 1] + (ldd + 1 + q) Cm[ldd, 0] c[k - 1][l, ldd - 1]
      - Sum[SeriesCoefficient[Tanh[x], {x, 0, j}] c[k - 1 - j][l, ldd], {j, 1, k - 1}])];
  c];

c0 = buildSeries2[0, q, 4];   (* lp = 0 *)
c2 = buildSeries2[2, q, 4];   (* lp = 2 *)

(* check against direct integration at small beta *)
bS = 1/20; etaS = ArcTanh[bS];
serVal[c_, l_, ldd_, qv_, ordE_] := Sum[(c[k][l, ldd] /. q -> qv) etaS^k, {k, 0, ordE}];
devTab = Table[{l, lp, ldd, qv,
    Abs[DD[qv, l, lp, ldd, 0, bS] - serVal[If[lp == 0, c0, c2], l, ldd, qv, 4]]},
   {lp, {0, 2}}, {l, 0, 3}, {ldd, 0, 3}, {qv, {1/2}}];
dev = Max[devTab[[All, All, All, All, 5]]];
Print["N8  ODE-generated Taylor series vs direct integration at beta=0.05 (O(eta^5)~3e-7 expected): max dev = ", ScientificForm[N[dev], 3]];

(* ---- N9: symbolic second-order tables (p = gamma beta expansion), diagonal elements ->
   check invariance under q -> -3-q (i.e. functions of Dhat_nu only) and print them. ---- *)
(* convert eta-series to p-series: p = sinh(eta) => eta = ArcSinh[p]; D as function of eta; compose *)
toP[c_, l_, ldd_, ordE_] := Module[{serEta, serP},
  serEta = Sum[c[k][l, ldd] et^k, {k, 0, ordE}];
  Normal@Series[serEta /. et -> ArcSinh[pp], {pp, 0, ordE}]];
Print["--- symbolic tables (m=0, q = -Ohat): ---"];
Do[Module[{expr = Collect[Simplify[toP[c0, l, ldd, 4]], pp]},
   Print["D_{", l, ",0,", ldd, "}(q) = ", expr]], {l, 0, 2}, {ldd, 0, 2}];
Do[Module[{expr = Collect[Simplify[toP[c2, l, ldd, 4]], pp]},
   Print["D_{", l, ",2,", ldd, "}(q) = ", expr]], {l, 0, 3}, {ldd, 0, 3}];
(* reflection invariance of diagonals *)
dev = Max@Table[Module[{e1 = toP[If[lp == 0, c0, c2], l, l, 4]},
    Max@Abs@CoefficientList[Expand[(e1 /. q -> -3 - q) - e1] /. pp -> zz, {zz, q}]],
   {lp, {0, 2}}, {l, 0, 3}];
report["N9  diagonal elements invariant under q -> -3-q (=> functions of Dhat_nu)", dev];

Print["--- doppler harness done ---"];
