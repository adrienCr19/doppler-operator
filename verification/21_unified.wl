(* ::Package:: *)
(* 21_unified.wl — the unified index-raising formulation.
   ONE object: the kernel column (F-bracket), a boost-operator element,
       F_{ab}(q) := ^{-1-q}K_{ab}(+beta)   (m = 0).
   ONE seed:   F_{00}(q) = J_{2+q}/(2p).
   TWO master recurrences (transposes of each other under the F-reflection): raise row a, or column b.
   The master identity: every Doppler element is bilinear in F,
       D_{l l' l''}(q) = ((-1)^{l'+l''}/gamma) F_{l l'}(q) F_{l' l''}(q-1).
   Symmetries: F-reflection F_{ab}(q) = (-1)^{a+b} F_{ba}(-4-q); J-evenness J_{-w}=J_w.
   Checks:
     V1  master identity vs direct Doppler definition
     V2a row-raise master recurrence for F
     V2b column-raise master recurrence for F
     V3  F-reflection
     V4  END-TO-END: build F up to indices 3 from the seed by the master recurrences ONLY,
         assemble D_313 (a new element), compare to the direct definition. *)

prec = 30;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{g = ga[b]}, 2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(g (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] /; l < 0 || lp < 0 := 0;
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
Cm0[l_] := If[l <= 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
JJ[w_, b_] := If[w == 0, 2 ArcTanh[b], With[{g = ga[b], p = ga[b] b}, ((g + p)^w - (g - p)^w)/w]];
report[name_, dev_, tol_: 10^-24] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (dev = ", ScientificForm[N[dev], 3], ")"];

bb = 3/10; qv = 1/2;

(* the one object, directly *)
Fdir[a_, b_, q_] := Kp[-1 - q, a, b, bb];
(* the physical Doppler operator, directly *)
Ddir[l_, lp_, ldd_, q_] := (1/ga[bb]) Kp[-1 - q, l, lp, bb] Km[-q, lp, ldd, bb];

(* ---- V1: master bilinear identity ---- *)
dev = Max@Table[Abs[Ddir[l, lp, ldd, qv] - ((-1)^(lp + ldd)/ga[bb]) Fdir[l, lp, qv] Fdir[lp, ldd, qv - 1]],
   {l, 0, 3}, {lp, 0, 3}, {ldd, 0, 3}];
report["V1  D_{l l' l''} = ((-1)^{l'+l''}/g) F_{l l'}(q) F_{l' l''}(q-1)", dev];

(* ---- V2a: row-raise master recurrence ---- *)
rowraise[a_, b_, q_] := (1/(Cm0[a] (ga[bb] bb))) (ga[bb] Fdir[a - 1, b, q] - Fdir[a - 1, b, q + 1]) - (Cm0[a - 1]/Cm0[a]) Fdir[a - 2, b, q];
dev = Max@Table[Abs[Fdir[a, b, qv] - rowraise[a, b, qv]], {a, 1, 3}, {b, 0, 3}];
report["V2a row-raise: F_{ab}(q) = (1/C_a p)(g F_{a-1,b}(q) - F_{a-1,b}(q+1)) - (C_{a-1}/C_a)F_{a-2,b}(q)", dev];

(* ---- V2b: column-raise master recurrence ---- *)
colraise[a_, b_, q_] := -(1/(Cm0[b] (ga[bb] bb))) (ga[bb] Fdir[a, b - 1, q] - Fdir[a, b - 1, q - 1]) - (Cm0[b - 1]/Cm0[b]) Fdir[a, b - 2, q];
dev = Max@Table[Abs[Fdir[a, b, qv] - colraise[a, b, qv]], {a, 0, 3}, {b, 1, 3}];
report["V2b col-raise: F_{ab}(q) = -(1/C_b p)(g F_{a,b-1}(q) - F_{a,b-1}(q-1)) - (C_{b-1}/C_b)F_{a,b-2}(q)", dev];

(* ---- V3: F-reflection ---- *)
dev = Max@Table[Abs[Fdir[a, b, qv] - (-1)^(a + b) Fdir[b, a, -4 - qv]], {a, 0, 3}, {b, 0, 3}];
report["V3  F-reflection: F_{ab}(q) = (-1)^{a+b} F_{ba}(-4-q)", dev];

(* ---- V4: END-TO-END generation from the seed alone ----
   Build Fgen from the seed F_{00}(q) = J_{2+q}/(2p) using ONLY the master recurrences (row then col),
   at the finite set of q-shifts required, then assemble D_313 and compare to Ddir. *)
seedF[q_] := JJ[2 + q, bb]/(2 ga[bb] bb);
ClearAll[Fgen];
Fgen[a_, b_, q_] /; a < 0 || b < 0 := 0;
Fgen[0, 0, q_] := seedF[q];
(* first fill column 0 (b=0) by row-raise, then raise columns by col-raise *)
Fgen[a_, 0, q_] := Fgen[a, 0, q] = (1/(Cm0[a] (ga[bb] bb))) (ga[bb] Fgen[a - 1, 0, q] - Fgen[a - 1, 0, q + 1]) - (Cm0[a - 1]/Cm0[a]) Fgen[a - 2, 0, q];
Fgen[a_, b_, q_] := Fgen[a, b, q] = -(1/(Cm0[b] (ga[bb] bb))) (ga[bb] Fgen[a, b - 1, q] - Fgen[a, b - 1, q - 1]) - (Cm0[b - 1]/Cm0[b]) Fgen[a, b - 2, q];
Dgen[l_, lp_, ldd_, q_] := ((-1)^(lp + ldd)/ga[bb]) Fgen[l, lp, q] Fgen[lp, ldd, q - 1];
Print["V4  END-TO-END from seed F00 only:"];
Do[report["     D_" <> ToString[t[[1]]] <> ToString[t[[2]]] <> ToString[t[[3]]] <> " generated vs direct",
    Abs[Dgen[t[[1]], t[[2]], t[[3]], qv] - Ddir[t[[1]], t[[2]], t[[3]], qv]]],
  {t, {{3, 1, 3}, {2, 1, 2}, {3, 0, 3}, {2, 2, 2}}}];

Print["--- unified formulation harness done ---"];
