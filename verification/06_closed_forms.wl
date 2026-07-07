(* ::Package:: *)
(* 06_closed_forms.wl — Verify NEW closed forms for higher-index Doppler operators,
   generated from the D000 base case (seed S(d) = ^dK_00(-beta)) via B3 + B10.
   Notation:  J_w = ((g+p)^w - (g-p)^w)/w   (odd part; J even in w),
              E_w = (g+p)^w + (g-p)^w       (even part; used after product collapse).
   Claims (m=0, s=0, on nu^q i.e. Ohat -> -q):
     K^d_10(-b)  = -(Sqrt[3]/(2p^2)) (g J_{2-d} - J_{1-d})     [from B10, one step]
     K^d_20(-b)  =  (Sqrt[5]/(4p^3)) ((3+2p^2) J_{1-d} - 6 g J_{2-d} + 3 J_{3-d})   [two steps; = stag240 9b]
     D101(q) = 3/(4 g p^4) (g J_{2+q} - J_{3+q}) (g J_{1+q} - J_{q})
     D020(q) = 5/(16 g p^6) ((3+2p^2) J_{2+q} - 6 g J_{1+q} + 3 J_{q})
                           * ((3+2p^2) J_{1+q} - 6 g J_{2+q} + 3 J_{3+q})
   Also verify the product-collapse rule J_{a+q}J_{b+q} = (E_{a+b+2q} - E_{a-b})/((a+q)(b+q)). *)

prec = 30;
nint[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, m_, x_] := Sqrt[(2 l + 1)/(4 Pi) Factorial[l - m]/Factorial[l + m]] LegendreP[l, m, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, m_, b_] := Module[{g = ga[b]},
   2 Pi nint[Ybar[l, m, mu] Ybar[lp, m, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, m_, b_] := Km[d, l, lp, m, b] = KmRaw[d, l, lp, m, b];
Kp[d_?NumericQ, l_, lp_, m_, b_] := Km[d, l, lp, m, -b];
DD[q_, l_, lp_, ldd_, b_] := (1/ga[b]) Kp[-1 - q, l, lp, 0, b] Km[-q, lp, ldd, 0, b];
report[name_, dev_, tol_: 10^-24] := Print[name, ": ", If[dev < tol, "PASS", "FAIL"], "  (max dev = ", ScientificForm[N[dev], 3], ")"];

(* J_w = ((g+p)^w-(g-p)^w)/w = 2 Sinh[w eta]/w  with eta = ArcTanh[b];  J_0 = 2 eta (limit) *)
JJ[w_, b_] := If[w == 0, 2 ArcTanh[b], With[{g = ga[b], p = ga[b] b}, ((g + p)^w - (g - p)^w)/w]];
EE[w_, b_] := With[{g = ga[b], p = ga[b] b}, ((g + p)^w + (g - p)^w)];

bs = {3/10, 7/10}; qsv = {-7/10, 1/2, 13/10, 2};

(* C1: kernel column closed forms *)
dev = Max@Table[With[{g = ga[b], p = ga[b] b},
    Abs[Km[d, 1, 0, 0, b] - (-(Sqrt[3]/(2 p^2)) (g JJ[1 - d, b] - JJ[2 - d, b]))]], {b, bs}, {d, {-2, -1/2, 1, 5/2}}];
report["F1  K^d_10(-b) = -Sqrt[3]/(2p^2) (g J_{1-d} - J_{2-d})", dev];
dev = Max@Table[With[{g = ga[b], p = ga[b] b},
    Abs[Km[d, 2, 0, 0, b] - (Sqrt[5]/(4 p^3)) ((3 + 2 p^2) JJ[1 - d, b] - 6 g JJ[2 - d, b] + 3 JJ[3 - d, b])]], {b, bs}, {d, {-2, -1/2, 1, 5/2}}];
report["F2  K^d_20(-b) closed form (matches stag240 9b)", dev];

(* C2: Doppler operator closed forms *)
dev = Max@Table[With[{g = ga[b], p = ga[b] b},
    Abs[DD[q, 1, 0, 1, b] - 3/(4 g p^4) (g JJ[2 + q, b] - JJ[3 + q, b]) (g JJ[1 + q, b] - JJ[q, b])]], {b, bs}, {q, qsv}];
report["F3  D101(q) = 3/(4 g p^4) (g J_{2+q}-J_{3+q})(g J_{1+q}-J_q)", dev];
dev = Max@Table[With[{g = ga[b], p = ga[b] b},
    Abs[DD[q, 0, 2, 0, b] - 5/(16 g p^6) ((3 + 2 p^2) JJ[2 + q, b] - 6 g JJ[1 + q, b] + 3 JJ[q, b]) *
       ((3 + 2 p^2) JJ[1 + q, b] - 6 g JJ[2 + q, b] + 3 JJ[3 + q, b])]], {b, bs}, {q, qsv}];
report["F4  D020(q) closed form", dev];

(* C3: base-case consistency: D000 = (1/(4 g p^2)) J_{2+q} J_{1+q} == CR26 (10a) *)
dev = Max@Table[With[{g = ga[b], p = ga[b] b},
    Abs[DD[q, 0, 0, 0, b] - 1/(4 g p^2) JJ[2 + q, b] JJ[1 + q, b]]], {b, bs}, {q, qsv}];
report["F5  D000 = J_{2+q}J_{1+q}/(4 g p^2)  (base case, = CR26 10a)", dev];

(* C4: product collapse J_{a+q}J_{b+q} = (E_{a+b+2q} - E_{a-b})/((a+q)(b+q)) *)
dev = Max@Table[Abs[JJ[a + q, b] JJ[bq + q, b] - (EE[a + bq + 2 q, b] - EE[a - bq, b])/((a + q) (bq + q))],
   {b, bs}, {q, qsv}, {a, {1, 2, 3}}, {bq, {0, 1, 2}}];
report["F6  J-product collapse to E functions", dev];

(* C5: reflection symmetry of the closed forms: q -> -3-q swaps the two brackets (diagonal elements) *)
dev = Max@Table[With[{g = ga[b], p = ga[b] b},
    Abs[3/(4 g p^4) (g JJ[2 + q, b] - JJ[3 + q, b]) (g JJ[1 + q, b] - JJ[q, b]) -
        (3/(4 g p^4) (g JJ[2 + qq, b] - JJ[3 + qq, b]) (g JJ[1 + qq, b] - JJ[qq, b]) /. qq -> -3 - q)]], {b, {3/10}}, {q, qsv}];
report["F7  closed-form D101 invariant under q -> -3-q", dev];

(* C6: series cross-check of D020 closed form against the ODE table:
   D020 = p^4 q(-12+5q+6q^2+q^3)/45 + O(p^6)  [Section 4 table, lp=2] *)
serTab[q_, p_] := p^4 q (-12 + 5 q + 6 q^2 + q^3)/45;
dev = Max@Table[With[{b = 1/20}, With[{g = ga[b], p = ga[b] b},
     Abs[5/(16 g p^6) ((3 + 2 p^2) JJ[2 + q, b] - 6 g JJ[1 + q, b] + 3 JJ[q, b]) *
          ((3 + 2 p^2) JJ[1 + q, b] - 6 g JJ[2 + q, b] + 3 JJ[3 + q, b]) - serTab[q, p]]]], {q, {1/2}}];
Print["F8  D020 closed form vs O(p^4) ODE series at beta=0.05 (expect ~p^6 ~ 1.6e-8): dev = ", ScientificForm[N[dev], 3]];

Print["--- closed-forms harness done ---"];

(* ---- Route 2 (novel N5/N6 recursions at the Doppler-operator level): intermediate steps ---- *)
(* lattice bases from seeds only: (j,k)D_000 = J_{2+j+q} J_{1+k+q} / (4 g p^2) *)
DL000[j_, k_, q_, b_] := With[{g = ga[b], p = ga[b] b}, JJ[2 + j + q, b] JJ[1 + k + q, b]/(4 g p^2)];
DLnum[j_, k_, q_, l_, lp_, ldd_, b_] := (1/ga[b]) Kp[-1 - j - q, l, lp, 0, b] Km[-k - q, lp, ldd, 0, b];
dev = Max@Table[Abs[DLnum[j, k, q, 0, 0, 0, 3/10] - DL000[j, k, q, 3/10]], {j, {0, 1}}, {k, {-1, 0}}, {q, {1/2, -7/10}}];
report["F9  lattice bases (j,k)D_000 = J_{2+j+q}J_{1+k+q}/(4gp^2)", dev];

(* N6 step: (j,k)D_001 = (Sqrt[3]/(4 g p^3)) J_{2+j+q} (g J_{1+k+q} - J_{k+q}) *)
DL001[j_, k_, q_, b_] := With[{g = ga[b], p = ga[b] b}, (Sqrt[3]/(4 g p^3)) JJ[2 + j + q, b] (g JJ[1 + k + q, b] - JJ[k + q, b])];
dev = Max@Table[Abs[DLnum[j, k, q, 0, 0, 1, 3/10] - DL001[j, k, q, 3/10]], {j, {0, 1}}, {k, {0}}, {q, {1/2, -7/10}}];
report["F10 N6 step: (j,0)D_001 closed form", dev];

(* N5 step assembles D101 = (1/(C1 p))[g D_001 - (1,0)D_001]  ->  must equal F3 form *)
dev = Max@Table[With[{g = ga[b], p = ga[b] b},
    Abs[(Sqrt[3]/p) (g DL001[0, 0, q, b] - DL001[1, 0, q, b]) -
        3/(4 g p^4) (g JJ[2 + q, b] - JJ[3 + q, b]) (g JJ[1 + q, b] - JJ[q, b])]], {b, {3/10, 7/10}}, {q, {1/2, -7/10}}];
report["F11 N5 assembly of D101 == kernel-route closed form (F3)", dev];
