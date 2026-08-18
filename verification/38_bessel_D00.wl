(* ::Package:: *)
(* 38_bessel_D00.wl — Rigorous Bessel identity for the thermally averaged D00 operator.

   Claim (proved in proofs/bessel-d00.pdf):

       <D000(q)>_theta  =  ( K_{2q+3}(1/theta) - K_1(1/theta) )
                           / ( 2 (q+1)(q+2) theta K_2(1/theta) ),

   where <.>_theta is the Maxwell-Juttner average <g> = Integral_0^inf p^2 fMB(p,theta) g(p) dp,
   fMB(p,theta) = Exp[-Sqrt[1+p^2]/theta]/(theta K_2(1/theta)), p = |p|/(m_e c) = Sinh[eta].

   This harness verifies, to deviation 0 where symbolic and to ~30 digits numerically,
   every non-trivial step of that proof. *)

Jw[w_, p_] := With[{eta = ArcSinh[p]}, If[w === 0, 2 eta, 2 Sinh[w eta]/w]];
D000c[q_, p_] := With[{g = Sqrt[1 + p^2]}, Jw[2 + q, p] Jw[1 + q, p]/(4 g p^2)];
fMB[p_?NumericQ, th_?NumericQ] := Exp[-Sqrt[1 + p^2]/th]/(th BesselK[2, 1/th]);

(* ---------- STEP 1: the integrand collapses -- measure x D000 has no rational prefactor ---------- *)
(* In rapidity eta (p = Sinh eta, dp = Cosh eta d eta), the measure is p^2 dp = Sinh^2 eta Cosh eta d eta.
   Multiplying by D000 = Jw[2+q] Jw[1+q]/(4 Cosh eta Sinh^2 eta) must cancel Cosh eta Sinh^2 eta exactly. *)
etaSym = Symbol["eta"];
measTimesD = Sinh[etaSym]^2 Cosh[etaSym] *
   (2 Sinh[(2 + q) etaSym]/(2 + q)) (2 Sinh[(1 + q) etaSym]/(1 + q)) /
   (4 Cosh[etaSym] Sinh[etaSym]^2);
collapsed = Sinh[(2 + q) etaSym] Sinh[(1 + q) etaSym]/((1 + q) (2 + q));
Print["S1  measure*D000 collapses to Sinh((2+q)e)Sinh((1+q)e)/((1+q)(2+q)): ",
  Simplify[measTimesD - collapsed] === 0];

(* ---------- STEP 2: product-to-sum ---------- *)
(* Sinh A Sinh B = 1/2 (Cosh(A+B) - Cosh(A-B)); A=(2+q)e, B=(1+q)e -> A+B=(3+2q)e, A-B=e *)
p2s = 1/2 (Cosh[(3 + 2 q) etaSym] - Cosh[etaSym]);
Print["S2  product-to-sum Sinh((2+q)e)Sinh((1+q)e) = 1/2(Cosh((3+2q)e)-Cosh(e)): ",
  Simplify[Sinh[(2 + q) etaSym] Sinh[(1 + q) etaSym] - p2s] === 0];

(* ---------- STEP 3: normalization Integral p^2 fMB dp = 1 (i.e. K_2(z)/z identity) ---------- *)
(* Integral_0^inf Sinh^2 e Cosh e Exp[-z Cosh e] d eta == K_2(z)/z, for z>0. Check numerically. *)
normInt[z_] := NIntegrate[Sinh[e]^2 Cosh[e] Exp[-z Cosh[e]], {e, 0, Infinity},
   WorkingPrecision -> 30, PrecisionGoal -> 20];
Print["S3  normalization Int Sinh^2 e Cosh e e^{-z Cosh e} de == K_2(z)/z: ",
  Table[ScientificForm[N[Abs[normInt[z] - BesselK[2, z]/z]], 3], {z, {5, 20, 50}}]];

(* ---------- STEP 4: Basset/Schlaefli integral representation K_nu(z)=Int e^{-z Cosh e}Cosh(nu e)de ---------- *)
(* verify for representative (possibly non-integer) orders used by the proof: nu = 3+2q *)
kRep[nu_, z_] := NIntegrate[Exp[-z Cosh[e]] Cosh[nu e], {e, 0, Infinity},
   WorkingPrecision -> 30, PrecisionGoal -> 20];
Print["S4  K_nu(z)=Int e^{-z Cosh e}Cosh(nu e)de for nu=1,2,3,4,5.5: ",
  Table[ScientificForm[N[Abs[kRep[nu, 12] - BesselK[nu, 12]]], 3], {nu, {1, 2, 3, 4, 11/2}}]];

(* ---------- STEP 5: the assembled closed form vs direct numerical average ---------- *)
avg[qv_, th_] := NIntegrate[p^2 fMB[p, th] D000c[qv, p], {p, 0, Infinity},
   WorkingPrecision -> 30, PrecisionGoal -> 18];
cf[qv_, th_] := (BesselK[2 qv + 3, 1/th] - BesselK[1, 1/th])/(2 (qv + 1) (qv + 2) th BesselK[2, 1/th]);
devs = Flatten@Table[Abs[avg[qv, th] - cf[qv, th]], {qv, {0, 1/2, 1, 3/2, -(1/2)}}, {th, {1/50, 1/10}}];
Print["S5  <D000(q)>_theta == closed form, max dev over {q,theta}: ",
  ScientificForm[N[Max[devs]], 3], "  ", If[Max[devs] < 10^-24, "PASS", "FAIL"]];

(* ---------- STEP 6: the q=0 consistency <D000(0)>=1 via the Bessel recurrence ---------- *)
(* K_3(z)-K_1(z) = (4/z) K_2(z) [from K_{v+1}-K_{v-1}=(2v/z)K_v, v=2], so cf(0)= (4/z K_2)/(2*2*(1/z)K_2)=1. *)
zz = 1/17;
Print["S6a K_3(z)-K_1(z) - (4/z)K_2(z) = 0: ",
  ScientificForm[Quiet@N[Abs[(BesselK[3, N[zz, 60]] - BesselK[1, N[zz, 60]]) - (4/zz) BesselK[2, N[zz, 60]]]], 3]];
Print["S6b closed form at q=0 equals exactly 1: ", N[cf[0, 1/20], 25]];
(* and D000(0,p)=1 pointwise, so the average is trivially 1 *)
Print["S6c D000(0,p) == 1 pointwise (symbolic): ",
  Simplify[D000c[0, p] - 1, Assumptions -> p > 0] === 0];

Print["--- 38_bessel_D00 done ---"];
