(* ::Package:: *)
(* 10_middle_index.wl — Middle-index (l') raising for the Doppler operator.
   Extended family (s=0):  X^{(j,k)}_{l; l1,l2; ldd}(q) = (1/g) K^{-1-j-q}_{l,l1}(+b) K^{-k-q}_{l2,ldd}(-b),
   physical Doppler operator = diagonal slice l1 = l2 = l'.
   Ingredient A (column recursion of the +b kernel; from B10 sandwiched in transposition B3):
     K^{a}_{l,lp}(+b) = -(1/(C_lp p)) [ g K^{a}_{l,lp-1}(+b) - K^{a+1}_{l,lp-1}(+b) ] - (C_{lp-1}/C_lp) K^{a}_{l,lp-2}(+b)
   Ingredient B (row recursion of the -b kernel; B10 as-is):
     K^{b}_{lp,ldd}(-b) = -(1/(C_lp p)) [ g K^{b}_{lp-1,ldd}(-b) - K^{b-1}_{lp-1,ldd}(-b) ] - (C_{lp-1}/C_lp) K^{b}_{lp-2,ldd}(-b)
   (M) product of A and B => diagonal step:
     D^{(j,k)}_{l,lp,ldd} =
        (1/(C_lp^2 p^2)) [ g^2 D^{(j,k)} - g D^{(j-1,k)} - g D^{(j,k+1)} + D^{(j-1,k+1)} ]_{l,lp-1,ldd}
      + (C_{lp-1}/(C_lp^2 p)) [ g X^{(j,k)}_{l;lp-1,lp-2;ldd} - X^{(j-1,k)}_{l;lp-1,lp-2;ldd} ]
      + (C_{lp-1}/(C_lp^2 p)) [ g X^{(j,k)}_{l;lp-2,lp-1;ldd} - X^{(j,k+1)}_{l;lp-2,lp-1;ldd} ]
      + (C_{lp-1}^2/C_lp^2) D^{(j,k)}_{l,lp-2,ldd} .
   All quantities evaluated on power laws nu^q by 30-digit quadrature. *)

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
report[name_, dev_, tol_: 10^-24] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

XX[j_, k_, q_, l_, l1_, l2_, ldd_, m_, b_] := (1/ga[b]) Kp[-1 - j - q, l, l1, m, b] Km[-k - q, l2, ldd, m, b];
DD[j_, k_, q_, l_, lp_, ldd_, m_, b_] := XX[j, k, q, l, lp, lp, ldd, m, b];

bb = 3/10; qs = {1/2, -7/10};

(* ---- MX1: Ingredient A, column recursion of the +beta kernel (general m) ---- *)
dev = Max@Table[With[{g = ga[bb], p = ga[bb] bb},
    Abs[Kp[a, l, lp, m, bb] - (-(1/(Cm[lp, m] p)) (g Kp[a, l, lp - 1, m, bb] - Kp[a + 1, l, lp - 1, m, bb]) -
        (Cm[lp - 1, m]/Cm[lp, m]) If[lp - 2 >= Abs[m], Kp[a, l, lp - 2, m, bb], 0])]],
   {a, {-3/2, 1/2, 2}}, {m, 0, 1}, {l, m, 2}, {lp, Max[Abs[m], 1] + 1, 3}];
report["M1  column recursion of K(+b)  [B10 transposed]", dev];

(* ---- M: the middle-index raising relation for the Doppler operator (m = 0 and m = 1) ---- *)
Mrel[j_, k_, q_, l_, lp_, ldd_, m_, b_] := With[{g = ga[b], p = ga[b] b, C1 = Cm[lp, m], C0 = Cm[lp - 1, m]},
   (1/(C1^2 p^2)) (g^2 DD[j, k, q, l, lp - 1, ldd, m, b] - g DD[j - 1, k, q, l, lp - 1, ldd, m, b] -
       g DD[j, k + 1, q, l, lp - 1, ldd, m, b] + DD[j - 1, k + 1, q, l, lp - 1, ldd, m, b]) +
    (C0/(C1^2 p)) (g XX[j, k, q, l, lp - 1, lp - 2, ldd, m, b] - XX[j - 1, k, q, l, lp - 1, lp - 2, ldd, m, b]) +
    (C0/(C1^2 p)) (g XX[j, k, q, l, lp - 2, lp - 1, ldd, m, b] - XX[j, k + 1, q, l, lp - 2, lp - 1, ldd, m, b]) +
    (C0^2/C1^2) If[lp - 2 >= Abs[m], DD[j, k, q, l, lp - 2, ldd, m, b], 0]];
dev = Max@Table[Abs[DD[0, 0, q, l, lp, ldd, m, bb] - Mrel[0, 0, q, l, lp, ldd, m, bb]],
   {q, qs}, {m, 0, 1}, {l, m, 2}, {lp, Max[Abs[m], 1] + 1, 3}, {ldd, m, 2}];
report["M2  middle-index raising (M) for the Doppler operator", dev];

(* ---- M3: two-step climb 0 -> 2 rebuilds D020 and matches the N8 closed form ---- *)
JJ[w_, b_] := If[w == 0, 2 ArcTanh[b], With[{g = ga[b], p = ga[b] b}, ((g + p)^w - (g - p)^w)/w]];
D020cf[q_, b_] := With[{g = ga[b], p = ga[b] b},
   5/(16 g p^6) ((3 + 2 p^2) JJ[2 + q, b] - 6 g JJ[1 + q, b] + 3 JJ[q, b]) *
                ((3 + 2 p^2) JJ[1 + q, b] - 6 g JJ[2 + q, b] + 3 JJ[3 + q, b])];
dev = Max@Table[Abs[Mrel[0, 0, q, 0, 2, 0, 0, bb] - D020cf[q, bb]], {q, qs}];
report["M3  (M) at lp=2 reproduces the D020 closed form", dev];

(* ---- M4: at lattice points (relation holds for all (j,k)) ---- *)
dev = Max@Table[Abs[DD[j, k, q, 1, 2, 1, 0, bb] - Mrel[j, k, q, 1, 2, 1, 0, bb]], {j, {0, 1}}, {k, {-1, 0}}, {q, {1/2}}];
report["M4  (M) at shifted lattice points (j,k)", dev];

Print["--- middle-index harness done ---"];

(* ---- M5: symbolic verification of the worked D020 rederivation (page N10, worked example).
   With A = g^2 J2 - 2 g J1 + J0,  B = g^2 J1 - 2 g J2 + J3  (lattice-block factors) the assembled (M) gives
   D020 = (5/(16 g p^6)) [ 9 A B - 3 p^2 (A J1 + J2 B) + p^4 J1 J2 ]
   which must equal the N8 closed form (5/(16 g p^6)) U V with
   U = (3+2p^2) J2 - 6 g J1 + 3 J0,  V = (3+2p^2) J1 - 6 g J2 + 3 J3,
   as a POLYNOMIAL identity in the free symbols J0..J3, using only g^2 = 1 + p^2. ---- *)
Clear[g, p, J0, J1, J2, J3];
A = g^2 J2 - 2 g J1 + J0; B = g^2 J1 - 2 g J2 + J3;
U = (3 + 2 p^2) J2 - 6 g J1 + 3 J0; V = (3 + 2 p^2) J1 - 6 g J2 + 3 J3;
lhs = 9 A B - 3 p^2 (A J1 + J2 B) + p^4 J1 J2;
diff = PolynomialReduce[Expand[lhs - U V], {g^2 - 1 - p^2}, {g, p, J0, J1, J2, J3}][[2]];
Print["M5  worked D020 rederivation: 9AB - 3p^2(A J1 + J2 B) + p^4 J1 J2 - UV  ==  ",
  If[diff === 0, "0 (PASS, exact polynomial identity mod g^2 = 1+p^2)", diff]];
