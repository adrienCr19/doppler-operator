(* ::Package:: *)
(* 12_selfadjoint.wl — N12: the Doppler operator is self-adjoint on photon phase space
   L^2(nu^2 dnu) x (sphere), and N1 is its coordinate expression.
   Key facts verified here on NON-power-law spectra:
     S1  Ohat^dagger = 3 - Ohat on L^2(nu^2 dnu):  <f, e^{a Ohat} g> = <e^{a(3-Ohat)} f, g>
     S2  diagonal element self-adjoint:   <f, D101 g> = <D101 f, g>
     S3  off-diagonal adjoint pair:       <f, D001 g> = <D100 f, g>
   Operator action via the integral representation  J_{w-Ohat} g(nu) = Int_{-eta}^{eta} e^{w u} g(e^u nu) du
   (from J_w = int_{-eta}^{eta} e^{wu} du and e^{c Ohat} g(nu) = g(e^{-c} nu)).
   Inner product: <f,g> = Int_0^inf f(nu) g(nu) nu^2 dnu.  *)

prec = 20;
b = 3/10; g0 = 1/Sqrt[1 - b^2]; p0 = g0 b; eta0 = ArcTanh[b];
f1[nu_] := Exp[-nu]; f2[nu_] := Exp[-2 nu];   (* non-power-law test spectra *)
report[name_, dev_, tol_: 10^-8] := Print[name, ": ", If[TrueQ[N[dev] < tol], "PASS", "FAIL"], "  (dev = ", ScientificForm[N[dev], 3], ")"];

(* ---- S1: Ohat adjoint ---- *)
a0 = 3/10;
lhs = NIntegrate[nu^2 f1[nu] f2[Exp[-a0] nu], {nu, 0, Infinity}, WorkingPrecision -> prec, PrecisionGoal -> 12];
rhs = NIntegrate[nu^2 Exp[3 a0] f1[Exp[a0] nu] f2[nu], {nu, 0, Infinity}, WorkingPrecision -> prec, PrecisionGoal -> 12];
report["S1  <f, e^{a Ohat} g> = <e^{a(3-Ohat)} f, g>", Abs[lhs - rhs], 10^-10];

(* ---- operator pairings via the 3D integral representation ----
   D101 = (3/(4 g p^4)) (g J_{2-O} - J_{3-O})(g J_{1-O} - J_{O})
   D001 = (Sqrt[3]/(4 g p^3)) J_{2-O} (g J_{1-O} - J_{O})
   D100 = (Sqrt[3]/(4 g p^3)) (g J_{2-O} - J_{3-O}) J_{1-O}
   <f, Op g> = pref * Int nu^2 f(nu) W1(u) W2(v) g(e^{u+v} nu) du dv dnu  *)
(* nu-integral done analytically for exponential test spectra:
   f = e^{-alpha nu}, g = e^{-beta nu}:  Int nu^2 f(nu) g(c nu) dnu = 2/(alpha + beta c)^3, c = e^{u+v} *)
pairE[W1_, W2_, pref_, alpha_, beta_] := pref NIntegrate[
   W1[u] W2[v] 2/(alpha + beta Exp[u + v])^3,
   {u, -eta0, eta0}, {v, -eta0, eta0},
   WorkingPrecision -> prec, PrecisionGoal -> 14, MaxRecursion -> 15];
W2m3[u_] := g0 Exp[2 u] - Exp[3 u];   (* g J_{2-O} - J_{3-O} *)
W1m0[v_] := g0 Exp[v] - 1;            (* g J_{1-O} - J_{O}   *)
Wj2[u_] := Exp[2 u];                  (* J_{2-O}             *)
Wj1[v_] := Exp[v];                    (* J_{1-O}             *)

(* S2: diagonal self-adjointness, f = e^{-nu}, g = e^{-2nu} *)
prefD101 = 3/(4 g0 p0^4);
lhs = pairE[W2m3, W1m0, prefD101, 1, 2];      (* <f, D101 g> *)
rhs = pairE[W2m3, W1m0, prefD101, 2, 1];      (* <D101 f, g> *)
report["S2  <f, D101 g> = <D101 f, g>  (non-power-law spectra)", Abs[lhs - rhs]/Abs[lhs], 10^-10];

(* S3: off-diagonal adjoint pair D001^dagger = D100 *)
prefD3 = Sqrt[3]/(4 g0 p0^3);
lhs = pairE[Wj2, W1m0, prefD3, 1, 2];         (* <f, D001 g> *)
rhs = pairE[W2m3, Wj1, prefD3, 2, 1];         (* <D100 f, g> *)
report["S3  <f, D001 g> = <D100 f, g>  (adjoint pair)", Abs[lhs - rhs]/Abs[lhs], 10^-10];

(* S4: sanity anchor — the same pairing evaluated on power laws reproduces the closed form:
   <nu^s e^{-nu}, D101 g> with g = nu^q ... instead simplest: apply operator rep to nu^q and
   compare with closed-form scalar (already F3); here check the rep itself on one power law *)
qv = 1/2;
JJ[w_] := If[w == 0, 2 eta0, ((g0 + p0)^w - (g0 - p0)^w)/w];
opOnPow = prefD101 NIntegrate[W2m3[u] W1m0[v] Exp[(u + v) qv], {u, -eta0, eta0}, {v, -eta0, eta0},
    WorkingPrecision -> prec, PrecisionGoal -> 12];
cf = prefD101 (g0 JJ[2 + qv] - JJ[3 + qv]) (g0 JJ[1 + qv] - JJ[qv]);
report["S4  integral rep on nu^q == closed form D101(q)", Abs[opOnPow - cf], 10^-12];

Print["--- self-adjointness harness done ---"];
