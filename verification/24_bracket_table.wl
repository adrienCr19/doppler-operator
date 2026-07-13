(* ::Package:: *)
(* 24_bracket_table.wl — verify the full bracket table and the closed-form assembly recipe:
   the EASIEST route to a closed form of any D_{hkl} from D000.
     Column bracket:  K_{ab}(q) = ^{-1-q}K_{ab}(+beta)   (built from seed by 2 recurrences)
     Assembly:        D_{hkl}(q) = ((-1)^{k+l}/gamma) K_{hk}(q) K_{kl}(q-1)
   Shorthand in the table: J_n means J_{n+q}. *)

prec = 30;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 15, MaxRecursion -> 20];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_?NumericQ, l_, lp_, b_] := Module[{gg = ga[b]}, 2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(gg (1 + b mu))^d]];
Km[d_?NumericQ, l_, lp_, b_] /; l < 0 || lp < 0 := 0;
Km[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_?NumericQ, l_, lp_, b_] := Km[d, l, lp, -b];
Cm0[l_] := If[l <= 0, 0, Sqrt[l^2/(4 l^2 - 1)]];
report[name_, dev_, tol_: 10^-22] := Print[name, ": ", If[TrueQ[N[dev, 30] < tol], "PASS", "FAIL"], "  (dev = ", ScientificForm[N[dev, 30], 3], ")"];

bb = 3/10; g = ga[bb]; p = ga[bb] bb;
J[w_] := If[w == 0, 2 ArcTanh[bb], ((g + p)^w - (g - p)^w)/w];   (* literal J at index w *)
Kdir[a_, b_, q_] := Kp[-1 - q, a, b, bb];                        (* the column, direct *)
Ddir[h_, k_, l_, q_] := (1/g) Kp[-1 - q, h, k, bb] Km[-q, k, l, bb];

(* ---------- the PUBLISHED bracket table (J_n = J_{n+q}) ---------- *)
Ftab[0, 0, q_] := J[2 + q]/(2 p);
Ftab[0, 1, q_] := Sqrt[3] (J[1 + q] - g J[2 + q])/(2 p^2);
Ftab[1, 0, q_] := Sqrt[3] (g J[2 + q] - J[3 + q])/(2 p^2);
Ftab[0, 2, q_] := Sqrt[5] (3 J[q] - 6 g J[1 + q] + (3 + 2 p^2) J[2 + q])/(4 p^3);
Ftab[2, 0, q_] := Sqrt[5] ((3 + 2 p^2) J[2 + q] - 6 g J[3 + q] + 3 J[4 + q])/(4 p^3);
Ftab[1, 1, q_] := 3 (g J[1 + q] - (2 + p^2) J[2 + q] + g J[3 + q])/(2 p^3);
Ftab[1, 2, q_] := Sqrt[15] (3 g J[q] - 3 (3 + 2 p^2) J[1 + q] + g (9 + 2 p^2) J[2 + q] - (3 + 2 p^2) J[3 + q])/(4 p^4);
Ftab[2, 1, q_] := Sqrt[15] ((3 + 2 p^2) J[1 + q] - g (9 + 2 p^2) J[2 + q] + 3 (3 + 2 p^2) J[3 + q] - 3 g J[4 + q])/(4 p^4);
Ftab[2, 2, q_] := 5 (3 (3 + 2 p^2) J[q] - 12 g (3 + p^2) J[1 + q] + 2 (27 + 24 p^2 + 2 p^4) J[2 + q] - 12 g (3 + p^2) J[3 + q] + 3 (3 + 2 p^2) J[4 + q])/(8 p^5);

qs = {1/2, -7/10, 13/10};
dev = Max@Table[Abs[Ftab[a, b, q] - Kdir[a, b, q]], {a, 0, 2}, {b, 0, 2}, {q, qs}];
report["T1  full bracket table (F_ab, a,b<=2) == kernel columns", dev];

(* ---------- the two recurrences build the table from the seed ---------- *)
ClearAll[Kgen];
Kgen[a_, b_, q_] /; a < 0 || b < 0 := 0;
Kgen[0, 0, q_] := J[2 + q]/(2 p);
Kgen[a_, 0, q_] := Kgen[a, 0, q] = (g Kgen[a - 1, 0, q] - Kgen[a - 1, 0, q + 1])/(Cm0[a] p) - (Cm0[a - 1]/Cm0[a]) Kgen[a - 2, 0, q];
Kgen[a_, b_, q_] := Kgen[a, b, q] = -(g Kgen[a, b - 1, q] - Kgen[a, b - 1, q - 1])/(Cm0[b] p) - (Cm0[b - 1]/Cm0[b]) Kgen[a, b - 2, q];
dev = Max@Table[Abs[Kgen[a, b, 1/2] - Kdir[a, b, 1/2]], {a, 0, 3}, {b, 0, 3}];
report["T2  brackets generated from seed by 2 recurrences == direct (a,b<=3)", dev];

(* ---------- assembly: D_hkl = ((-1)^{k+l}/g) K_hk(q) K_kl(q-1) ---------- *)
Dfrom[h_, k_, l_, q_] := ((-1)^(k + l)/g) Kgen[h, k, q] Kgen[k, l, q - 1];
dev = Max@Table[Abs[Dfrom[h, k, l, 1/2] - Ddir[h, k, l, 1/2]], {h, 0, 3}, {k, 0, 3}, {l, 0, 3}];
report["T3  assembly D_hkl == direct (all h,k,l<=3)", dev];

(* ---------- worked CLOSED FORM (symbolic in q): D_313 and D_121 ---------- *)
Clear[qq]; Js[w_] := ((gg + pp)^w - (gg - pp)^w)/w;   (* symbolic J with gg=gamma, pp=p *)
KgenS[a_, b_] /; a < 0 || b < 0 := 0;
KgenS[0, 0] := Js[2 + qq]/(2 pp);
(* build symbolic brackets at the needed q-shifts by direct table (already verified) using symbolic J *)
FtabS[0, 0, s_] := Js[2 + qq + s]/(2 pp);
FtabS[1, 2, s_] := Sqrt[15] (3 gg Js[qq + s] - 3 (3 + 2 pp^2) Js[1 + qq + s] + gg (9 + 2 pp^2) Js[2 + qq + s] - (3 + 2 pp^2) Js[3 + qq + s])/(4 pp^4);
FtabS[2, 1, s_] := Sqrt[15] ((3 + 2 pp^2) Js[1 + qq + s] - gg (9 + 2 pp^2) Js[2 + qq + s] + 3 (3 + 2 pp^2) Js[3 + qq + s] - 3 gg Js[4 + qq + s])/(4 pp^4);
D121S = Simplify[((-1)^(2 + 1)/gg) FtabS[1, 2, 0] FtabS[2, 1, -1]];
Print["T4  D121 closed form (symbolic, via brackets):"];
Print["    ", InputForm[D121S]];
(* numeric cross-check of the symbolic form *)
dev = Abs[(D121S /. {gg -> g, pp -> p, qq -> 1/2}) - Ddir[1, 2, 1, 1/2]];
report["T4b D121 symbolic-bracket form == direct at q=1/2", dev];

Print["--- bracket table harness done ---"];
