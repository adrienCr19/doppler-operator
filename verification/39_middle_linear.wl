(* ::Package:: *)
(* 39_middle_linear.wl — The middle-index recurrence as an EXPLICIT LINEAR relation with
   ELEMENTARY (J-function) coefficients, for GENERAL l'.  Improves on P3 (proofs/recurrences.html),
   whose coefficients were written as ratios of Doppler operators: here every coefficient is a
   ratio of the finite axis brackets J_a(q), i.e. a closed form in J_w, g, p only.  No kernels,
   no extended family, no weight lattice appear in any coefficient.

   Result (border l''=0, any outer h, l'>=1):
     D_{h,l',0}(q) = a_l'(q) D_{h,l'-1,0}(q) + b_l'(q) D_{h,l'-1,0}(q-1) + c_l'(q) D_{h,l'-2,0}(q),
       a_l'(q) =  g(2l'+1)/(l' p)          J_{l'}(q-1)/J_{l'-1}(q-1),
       b_l'(q) = -(2l'+1)/(l' p)            J_{l'}(q-1)/J_{l'-1}(q-2),
       c_l'(q) = -(l'-1)(2l'+1)/(l'(2l'-3)) J_{l'}(q-1)/J_{l'-2}(q-1)   (c_1 = 0),
   where the axis bracket
     J_a(q) = gamma^{1+q} Int_{-1}^{1} P_a(mu)(1 - beta mu)^{1+q} dmu
            = (-1)^a Sum_{n=0}^a [mu^n]P_a(mu) * T_n(1+q),
     T_n(s) = p^{-(n+1)} Sum_{k=0}^n Binomial[n,k] (-g)^k J_{s+n+1-k},
   and the axis elements themselves are D_{a00}(q) = Sqrt(2a+1)/(4 g p) * J_{1+q} * J_a(q).
   On nu^q, m=0. *)

prec = 40;
nintg[f_] := NIntegrate[f, {mu, -1, 1}, WorkingPrecision -> prec, PrecisionGoal -> 22, MaxRecursion -> 26];
Ybar[l_, x_] := Sqrt[(2 l + 1)/(4 Pi)] LegendreP[l, x];
ga[b_] := 1/Sqrt[1 - b^2];
KmRaw[d_, l_, lp_, b_] := Module[{gg = ga[b]},
   2 Pi nintg[Ybar[l, mu] Ybar[lp, (mu + b)/(1 + b mu)]/(gg (1 + b mu))^d]];
Km[d_, l_, lp_, b_] /; l < 0 || lp < 0 := 0;
Km[d_, l_, lp_, b_] := Km[d, l, lp, b] = KmRaw[d, l, lp, b];
Kp[d_, l_, lp_, b_] := Km[d, l, lp, -b];
report[name_, dev_, tol_: 10^-12] := Print[name, ": ",
   If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (dev = ", ToString[N[dev], InputForm], ")"];

bb = 42/100; g = ga[bb]; p = ga[bb] bb;
JJ[w_] := If[w == 0, 2 ArcTanh[bb], ((g + p)^w - (g - p)^w)/w];
(* independent reference, built from kernels — used ONLY to check *)
Ddir[l_, lp_, ldd_, q_] := (1/g) Kp[-1 - q, l, lp, bb] Km[-q, lp, ldd, bb];

(* ---- (1) moment closed form: T_n(s) = g^s Int mu^n (1+b mu)^s dmu ---- *)
Tclosed[n_, s_] := p^(-(n + 1)) Sum[Binomial[n, k] (-g)^k JJ[s + n + 1 - k], {k, 0, n}];
Tnum[n_, s_] := g^s nintg[mu^n (1 + bb mu)^s];
report["T1  moment closed form T_n(s) vs integral (n<=5, s in {1.37,-0.6})",
   Max@Table[Abs[Tclosed[n, s] - Tnum[n, s]], {n, 0, 5}, {s, {137/100, -6/10}}]];

(* ---- (2) axis bracket J_a(q) and axis elements D_{a00} ---- *)
Jbr[a_, q_] := (-1)^a Sum[Coefficient[LegendreP[a, x], x, n] Tclosed[n, 1 + q], {n, 0, a}];
Da00[a_, q_] := (Sqrt[2 a + 1]/(4 g p)) JJ[1 + q] Jbr[a, q];
report["T2  axis element D_{a00} = Sqrt(2a+1)/(4gp) J_{1+q} J_a(q) vs kernel (a<=5)",
   Max@Table[Abs[Da00[a, q] - Ddir[a, 0, 0, q]], {a, 0, 5}, {q, {37/100, -(13/10)}}]];

(* ---- (3) the GENERAL-l' explicit LINEAR middle recurrence ---- *)
aG[lp_, q_] := (g (2 lp + 1)/(lp p)) Jbr[lp, q - 1]/Jbr[lp - 1, q - 1];
bG[lp_, q_] := -((2 lp + 1)/(lp p)) Jbr[lp, q - 1]/Jbr[lp - 1, q - 2];
cG[lp_, q_] := If[lp >= 2, -((lp - 1) (2 lp + 1)/(lp (2 lp - 3))) Jbr[lp, q - 1]/Jbr[lp - 2, q - 1], 0];
(* generate the whole border tower from the axis row using ONLY these J-explicit scalars *)
Dg[h_, 0, q_] := Da00[h, q];
Dg[h_, lp_, q_] /; lp < 0 := 0;
Dg[h_, lp_, q_] := Dg[h, lp, q] = aG[lp, q] Dg[h, lp - 1, q] + bG[lp, q] Dg[h, lp - 1, q - 1] + cG[lp, q] Dg[h, lp - 2, q];
report["T3  general-l' J-explicit LINEAR middle recurrence -> D_{h,l',0} (h<=3, l'<=5)",
   Max@Table[Abs[Dg[h, lp, q] - Ddir[h, lp, 0, q]], {h, 0, 3}, {lp, 1, 5}, {q, {37/100, -(13/10)}}]];

(* ---- (4) the l'=1 showcase: a clean two-term step ---- *)
(* D_{h,1,0}(q) = (3/p^2)(g J_{1+q} - J_{2+q})[ (g/J_{1+q}) D_{h00}(q) - (1/J_q) D_{h00}(q-1) ] *)
mid1[h_, q_] := (3/p^2) (g JJ[1 + q] - JJ[2 + q]) ((g/JJ[1 + q]) Da00[h, q] - (1/JJ[q]) Da00[h, q - 1]);
report["T4  l'=1 showcase two-term step (h<=5)",
   Max@Table[Abs[mid1[h, q] - Ddir[h, 1, 0, q]], {h, 0, 5}, {q, {37/100, -(13/10)}}]];

(* ---- (5) the sqrt/C_l prefactors collapse to the stated rationals (symbolic) ---- *)
CCs[l_] := Sqrt[l^2/(4 l^2 - 1)];
sA = Simplify[Table[(1/CCs[l]) Sqrt[(2 l + 1)/(2 l - 1)] - (2 l + 1)/l, {l, 1, 6}]];
sC = Simplify[Table[(CCs[l - 1]/CCs[l]) Sqrt[(2 l + 1)/(2 l - 3)] - ((l - 1) (2 l + 1))/(l (2 l - 3)), {l, 2, 6}]];
Print["T5  prefactor collapse (a,c): ", If[Union[sA] === {0} && Union[sC] === {0}, "PASS (all 0)", {sA, sC}]];

(* ---- (6) the ELEMENTARY <-> AXIS boundary: off the axis the outer raise is NOT elementary ---- *)
(* using the axis coefficient J_{1+q}/J_{2+q} to raise the outer index OFF the axis fails by O(1) *)
outAxisGuess[h_, lp_, ldd_, q_] := (1/CCs[h]) ((g/p) Ddir[h - 1, lp, ldd, q]
   - (1/p) (JJ[1 + q]/JJ[2 + q]) Ddir[h - 1, lp, ldd, q + 1] - CCs[h - 1] Ddir[h - 2, lp, ldd, q]);
errOff = Max@Table[Abs[outAxisGuess[h, 1, 1, q] - Ddir[h, 1, 1, q]], {h, 2, 4}, {q, {37/100}}];
Print["T6  outer raise with AXIS coeff OFF-axis (l'=l''=1): err = ",
   ToString[N[errOff], InputForm], "  => genuinely NON-elementary off the axis (P2 rank-one fill is the tool there)"];

(* ---- (7) the spectral shift q->q-1 is ESSENTIAL: no fixed-q middle recurrence exists ---- *)
(* Is D_{.,L,0}(q) in the h-span of the lower FIXED-q levels {D_{.,j,0}(q): j<L}?  (least squares) *)
hs = Range[0, 9]; qz = 37/100;
fixedqResid[L_] := Module[{A, b, c}, A = Table[Ddir[h, j, 0, qz], {h, hs}, {j, 0, L - 1}];
   b = Table[Ddir[h, L, 0, qz], {h, hs}]; c = LeastSquares[N[A, 30], N[b, 30]];
   Norm[N[A, 30].c - N[b, 30]]/Norm[N[b, 30]]];
(* ...and once the q-1 levels are adjoined (the P3 shift), it becomes exactly dependent *)
shiftResid[L_] := Module[{A, b, c}, A = Table[Flatten@Table[{Ddir[h, j, 0, qz], Ddir[h, j, 0, qz - 1]}, {j, 0, L - 1}], {h, hs}];
   b = Table[Ddir[h, L, 0, qz], {h, hs}]; c = LeastSquares[N[A, 30], N[b, 30]];
   Norm[N[A, 30].c - N[b, 30]]/Norm[N[b, 30]]];
fq = Max@Table[fixedqResid[L], {L, 2, 4}]; sh = Max@Table[shiftResid[L], {L, 2, 4}];
Print["T7  fixed-q middle dependence: worst residual = ", ToString[N[fq, 4]],
   "  (O(1) => NO fixed-q linear middle recurrence)"];
Print["    with the q->q-1 shift adjoined: worst residual = ", ToString[N[sh], InputForm],
   "  (=> exactly the P3 relation; the spectral shift is essential)"];

Print["--- middle-linear harness done ---"];
