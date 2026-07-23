(* ::Package:: *)
(* 37_identity_reductions.wl — SYMBOLIC (algebraic) proof that each Doppler-operator identity, once the
   operators are written out in terms of the aberration kernels, reduces to a kernel lemma or a tautology.
   Kernels are abstract commuting symbols kp[d,l,lp] = ^d K_{l,lp}(+beta),  km[d,l,lp] = ^d K_{l,lp}(-beta).
   Doppler operator (power-law reduction):  D(l,lp,ldd,q) = (1/g) kp[-1-q,l,lp] km[-q,lp,ldd].
   Kernel lemmas used as rewrite rules (each proved from the integral in harness 36):
     parity        km[d,l,lp] = (-1)^(l+lp) kp[d,l,lp]
     seed          kp[d,0,0]  = J[1-d]/(2p)   (and km via parity)
     transposition kp[d,l,lp] = km[2-d,lp,l]
     master (KR+)  kp[d,l,lp] = 1/(CC[l] p)(g kp[d,l-1,lp]-kp[d-1,l-1,lp]) - CC[l-1]/CC[l] kp[d,l-2,lp]
     column  (A)   kp[a,l,lp] = -1/(CC[lp] p)(g kp[a,l,lp-1]-kp[a+1,l,lp-1]) - CC[lp-1]/CC[lp] kp[a,l,lp-2]
   CC[.] and J[.] are inert heads. *)

ClearAll[kp, km, Dop, CC, J];
Dop[l_, lp_, ldd_, q_] := (1/g) kp[-1 - q, l, lp] km[-q, lp, ldd];
report[name_, expr_] := Module[{s = Simplify[expr]}, Print[name, ": ", If[s === 0, "0  -> VALID (LHS = RHS algebraically)", InputForm[s]]]];

(* ---- rewrite rules (kernel lemmas) ---- *)
parity = km[d_, l_, lp_] :> (-1)^(l + lp) kp[d, l, lp];
seed   = kp[d_, 0, 0] :> J[1 - d]/(2 p);
trKp   = kp[dd_, a_, b_] :> km[2 - dd, b, a];

(* ================= N1 reflection:  D_{h lp l}(q) = D_{l lp h}(-3-q) ================= *)
(* expand both sides, push every kp to km by transposition; the two become equal (km commute). *)
n1 = Dop[h, lp, l, q] - Dop[l, lp, h, -3 - q];
report["N1  reflection (via transposition)", n1 /. trKp];

(* ================= P2 rank-one:  D_{h lp l} = D_{h lp 0} D_{0 lp l} / D_{0 lp 0} ================= *)
(* pure cancellation of the shared kp[-1-q,0,lp] and km[-q,lp,0]; no lemma needed. *)
p2 = Dop[h, lp, l, q] - Dop[h, lp, 0, q] Dop[0, lp, l, q]/Dop[0, lp, 0, q];
report["P2  rank-one fill (cancellation)", Together[p2]];

(* ================= P0 outer axis ladder ================= *)
(* claim in D; expand; seed reduces the factor-2 monopole; result is the master recursion (KR+). *)
p0 = ( CC[a] Dop[a, 0, 0, q]
   - ((g/p) Dop[a - 1, 0, 0, q] - (1/p) (J[1 + q]/J[2 + q]) Dop[a - 1, 0, 0, q + 1] - CC[a - 1] Dop[a - 2, 0, 0, q]) );
p0seed = p0 /. km[d_, 0, 0] :> J[1 - d]/(2 p);      (* factor-2 monopoles -> J's *)
masterA = kp[-1 - q, a, 0] :> (1/(CC[a] p)) (g kp[-1 - q, a - 1, 0] - kp[-2 - q, a - 1, 0]) - (CC[a - 1]/CC[a]) kp[-1 - q, a - 2, 0];
report["P0  outer axis ladder (-> master recursion KR+)", p0seed /. masterA];

(* ================= P3 middle-index step ================= *)
(* claim in D; expand every D (incl. the axis coefficients); parity+seed collapse the factor-2 monopoles;
   Together cancels the shared axis kernels; substituting the column recursion (A) for the raised kernel
   kp[-1-q,h,lp+1] gives 0.  (lp,h declared integer so (-1)^(lp+1) = -(-1)^lp etc.) *)
p3 = ( CC[lp + 1] Dop[h, lp + 1, 0, q]
   - ( (Dop[lp + 1, 0, 0, q - 1]/p) (g Dop[h, lp, 0, q]/Dop[lp, 0, 0, q - 1]
        - (J[q - 1]/J[q]) Dop[h, lp, 0, q - 1]/Dop[lp, 0, 0, q - 2])
      - CC[lp] (Dop[lp + 1, 0, 0, q - 1]/Dop[lp - 1, 0, 0, q - 1]) Dop[h, lp - 1, 0, q] ) );
colA = kp[-1 - q, h, lp + 1] :>
   -(1/(CC[lp + 1] p)) (g kp[-1 - q, h, lp] - kp[-q, h, lp]) - (CC[lp]/CC[lp + 1]) kp[-1 - q, h, lp - 1];
p3k = Together[ p3 /. parity /. seed ];          (* all kernels -> kp; monopoles -> J; cancel *)
p3r = p3k /. colA;                               (* raised kernel -> column recursion (A) *)
report["P3  middle step (-> column recursion A)", p3r];

Print["--- identity-reduction harness done ---"];
