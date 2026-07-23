(* ::Package:: *)
(* 35_rederive_closed_forms.wl — rederive the closed forms of D020 and D121 using ONLY the N25
   pure-Doppler relations (P0 axis ladder, P3 middle border step, N1 reflection, P2 rank-one fill),
   starting from the seed D000.  Two checks each:
     (S) SYMBOLIC exact identity vs the published N8/N16 closed forms, proved by mapping every
         J_w -> (x^{w-q} y - x^{q-w}/y)/w  (x = e^eta, y = e^{q eta}) and Together-ing to 0.
     (N) NUMERIC vs direct kernel quadrature at beta in {0.3,0.7}, q in {1/2,-7/10}. *)

(* ---------- symbolic layer: J[w], g, p are inert; q, x, y symbolic ---------- *)
CC[l_] := Sqrt[l^2/(4 l^2 - 1)];                     (* C_l *)
J[w_] := JJ[w];                                       (* inert head *)
D000[q_] := JJ[2 + q] JJ[1 + q]/(4 g p^2);            (* seed *)

(* P0 axis ladders (pure D) *)
ClearAll[Da00, D00b];
Da00[0, q_] := D000[q];  Da00[a_, q_] /; a < 0 := 0;
Da00[a_, q_] := Da00[a, q] = Together[
   (1/CC[a]) ((g/p) Da00[a - 1, q] - (1/p) (JJ[1 + q]/JJ[2 + q]) Da00[a - 1, q + 1] - CC[a - 1] Da00[a - 2, q])];
D00b[0, q_] := D000[q];  D00b[b_, q_] /; b < 0 := 0;
D00b[b_, q_] := D00b[b, q] = Together[
   (1/CC[b]) ((g/p) D00b[b - 1, q] - (1/p) (JJ[2 + q]/JJ[1 + q]) D00b[b - 1, q - 1] - CC[b - 1] D00b[b - 2, q])];
Dax[a_, q_] := If[a == 0, D000[q], Da00[a, q]];

(* P3 middle-index self-fed border tower (pure D): builds D_{h,l',0} *)
ClearAll[Dgen];
Dgen[h_, 0, q_] := Dax[h, q];
Dgen[h_, lp_, q_] /; lp < 0 := 0;
Dgen[h_, lp_, q_] := Dgen[h, lp, q] = Together[
   (1/CC[lp]) (
      (Dax[lp, q - 1]/p) (
          g  Dgen[h, lp - 1, q]/Dax[lp - 1, q - 1]
        - (JJ[q - 1]/JJ[q]) Dgen[h, lp - 1, q - 1]/Dax[lp - 1, q - 2] )
      - If[lp >= 2, CC[lp - 1] (Dax[lp, q - 1]/Dax[lp - 2, q - 1]) Dgen[h, lp - 2, q], 0] )];

(* ---------- the two targets, purely from N25 ---------- *)
D020tower[q_] := Dgen[0, 2, q];                                   (* border, l''=0 *)
(* interior via N1 (2nd border) + P2 fill *)
D121tower[q_] := Dgen[1, 2, q] Dgen[1, 2, -3 - q]/Dgen[0, 2, q];

(* ---------- published closed forms ---------- *)
Apub[q_] := (3 + 2 p^2) JJ[2 + q] - 6 g JJ[1 + q] + 3 JJ[q];
Bpub[q_] := (3 + 2 p^2) JJ[1 + q] - 6 g JJ[2 + q] + 3 JJ[3 + q];
D020pub[q_] := (5/(16 g p^6)) Apub[q] Bpub[q];
Wpub[q_] := 3 g JJ[q] - 3 (3 + 2 p^2) JJ[1 + q] + g (9 + 2 p^2) JJ[2 + q] - (3 + 2 p^2) JJ[3 + q];
D121pub[q_] := (15/(16 g p^8)) Wpub[q] Wpub[-3 - q];

(* ---------- SYMBOLIC identity check ---------- *)
(* faithful map to two independent exponentials  X = e^eta,  T = e^{q eta}:
   every J-argument is integer +/- q, so  J[w] = (X^{n} T^{c} - X^{-n} T^{-c})/w
   with n = (w at q->0), c = Coefficient[w,q] in {+1,-1}.  X,T,q treated free. *)
jToXT[w_] := (X^(w /. q -> 0) T^Coefficient[w, q] - X^(-(w /. q -> 0)) T^(-Coefficient[w, q]))/w;
toExp = { JJ[w_] :> jToXT[w], g -> (X + 1/X)/2, p -> (X - 1/X)/2 };
zeroQ[e_] := Together[PowerExpand[e /. toExp]] // Numerator // Expand;

brief[e_] := If[e === 0, "0  ->  IDENTICAL (PASS)", Row[{"NONZERO (FAIL): ", Length[e], " terms remain"}]];
s020 = zeroQ[D020tower[q] - D020pub[q]];
Print["S-D020  tower(N25) - published(N8 eq.1.3)  ->  ", brief[s020]];
s121 = zeroQ[D121tower[q] - D121pub[q]];
Print["S-D121  tower(N25) - published(N16 box)    ->  ", brief[s121]];

(* also: does the tower-derived D020 factor into A(q)B(q)?  print the simplified prefactor identity *)
Print["S-D020f (16 g p^6 / 5) D020tower / (A B)   ->  ",
   Together[PowerExpand[((16 g p^6/5) D020tower[q]/(Apub[q] Bpub[q])) /. toExp]]];

(* ---------- NUMERIC cross-check vs direct kernel quadrature ---------- *)
prec = 40;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 20, MaxRecursion -> 24];
Ybar[l_, xx_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, xx];
gaN[b_] := 1/Sqrt[1 - b^2];
KmN[d_?NumericQ, l_, lp_, b_] := KmN[d, l, lp, b] =
   2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(gaN[b] (1 + b mu))^d];
KpN[d_?NumericQ, l_, lp_, b_] := KmN[d, l, lp, -b];
DdirN[l_, lp_, ldd_, q_, b_] := (1/gaN[b]) KpN[-1 - q, l, lp, b] KmN[-q, lp, ldd, b];

numChk[b_, q0_] := Module[{gg = gaN[b], pp = gaN[b] b, rule, t020, t121},
   rule = { g -> gg, p -> pp, JJ[w_] :> ((gg + pp)^w - (gg - pp)^w)/w };
   t020 = D020tower[q0] /. rule;  t121 = D121tower[q0] /. rule;
   { Abs[t020 - DdirN[0, 2, 0, q0, b]], Abs[t121 - DdirN[1, 2, 1, q0, b]] }];

tol = 10^-25;
Do[ With[{r = numChk[b, q0]},
     Print["N  beta=", N[b], " q=", N[q0],
        "  D020 ", If[N[r[[1]]] < tol, "PASS", "FAIL"], " (", ToString[N[r[[1]], 2], InputForm], ")",
        "  D121 ", If[N[r[[2]]] < tol, "PASS", "FAIL"], " (", ToString[N[r[[2]], 2], InputForm], ")"] ],
   {b, {3/10, 7/10}}, {q0, {1/2, -7/10}}];

Print["--- rederivation harness done ---"];
