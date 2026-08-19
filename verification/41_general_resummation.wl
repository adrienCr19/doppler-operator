(* ::Package:: *)
(* 41_general_resummation.wl — CLOSED-FORM SOLUTION of the Doppler recurrence ladder.
   Instead of stepping the ladder element by element, the middle-index recurrence (P3') is
   SUMMED IN CLOSED FORM, giving a general formula for arbitrary (l,l',l'') directly in terms
   of the AXIS family D_{a00} alone (which itself comes from the seed D000 by P0).

   Normalised axis element (removes the monopole seed factor):
       Dhat_a(q) := D_{a00}(q) / J_{q+1}
   RESUMMATION (theorem):
       D_{h,l',0}(q) = A_{l'}(q) * Sum_{j=0}^{l'} c^{(l')}_j * Dhat_h(q-j),
   where the c^{(l')}_j are ELEMENTARY in gamma,p (no J at all) and are independent of the outer
   index h, and A_{l'}(q) is an h-independent prefactor.  The c's turn out to be exactly the
   coefficients of the axis bracket, e.g. c^(2) ~ {3+2p^2, -6 gamma, 3} — the same coefficients
   that appear in the published D020 closed form.
   MASTER FORMULA (with P2 rank-one fill + N1 reflection):
       D_{l l' l''}(q) = [Sum_j W_j(q)   D_{l00}(q-j)]
                       * [Sum_k W_k(-3-q) D_{l''00}(-3-q-k)] / [Sum_m W_m(q) D_{000}(q-m)]
   No kernels anywhere; on nu^q, m=0.  (Operator form: q -> -Ohat.) *)

g=(X+1/X)/2; p=(X-1/X)/2;                      (* X = e^eta, internal only *)
sh[e_,d_]:=e/.Jf[m_]:>Jf[m+d];                 (* q -> q+d *)
refl[e_]:=e/.Jf[m_]:>Jf[3-m];                  (* N1: q -> -3-q *)
CC[l_]:=Sqrt[l^2/(4 l^2-1)]; red[e_]:=Cancel[Together[e]];
D000=Jf[2]Jf[1]/(4 g p^2); Da[0]=D000; Da[a_]/;a<0:=0;
Da[a_]:=Da[a]=red[(1/CC[a])((g/p)Da[a-1]-(1/p)(Jf[1]/Jf[2])sh[Da[a-1],1]-CC[a-1]Da[a-2])];
Tn[n_]:=(1/p^(n+1))Sum[Binomial[n,k](-g)^k Jf[n+2-k],{k,0,n}];
Jbr[a_]:=(-1)^a Sum[Coefficient[LegendreP[a,x],x,n]Tn[n],{n,0,a}];
aC[lp_]:=(g(2 lp+1)/(lp p)) sh[Jbr[lp],-1]/sh[Jbr[lp-1],-1];
bC[lp_]:=-((2 lp+1)/(lp p)) sh[Jbr[lp],-1]/sh[Jbr[lp-1],-2];
cC[lp_]:=If[lp>=2,-((lp-1)(2 lp+1)/(lp(2 lp-3))) sh[Jbr[lp],-1]/sh[Jbr[lp-2],-1],0];
Db[h_,0]:=Da[h]; Db[h_,lp_]/;lp<0:=0;
Db[h_,lp_]:=Db[h,lp]=red[aC[lp]Db[h,lp-1]+bC[lp]sh[Db[h,lp-1],-1]+cC[lp]Db[h,lp-2]];
Dfull[l_,lp_,ldd_]:=red[Db[l,lp] refl[Db[ldd,lp]]/Db[0,lp]];

(* ---- (1) unroll P3' into shift weights W^{(l')}_j ---- *)
W[0]={1}; W[lp_]/;lp<0:={};
W[lp_]:=W[lp]=Module[{Wm1=W[lp-1],Wm2=If[lp>=2,W[lp-2],{}],out},
  out=ConstantArray[0,lp+1];
  Do[out[[j]]=out[[j]]+aC[lp] Wm1[[j]],{j,1,Length[Wm1]}];
  Do[out[[j+1]]=out[[j+1]]+bC[lp] sh[Wm1[[j]],-1],{j,1,Length[Wm1]}];
  If[lp>=2,Do[out[[j]]=out[[j]]+cC[lp] Wm2[[j]],{j,1,Length[Wm2]}]];
  red/@out];
Bord[h_,lp_]:=Sum[W[lp][[j]] sh[Da[h],-(j-1)],{j,1,Length[W[lp]]}];
Dmaster[l_,lp_,ldd_]:=red[Bord[l,lp] refl[Bord[ldd,lp]]/Bord[0,lp]];

(* ---- numeric reference: the aberration kernel (CHECK ONLY) ---- *)
bb=42/100; gN=1/Sqrt[1-bb^2]; etaN=ArcSinh[gN bb]; qN=37/100;
ev[e_]:=N[e/.X->Exp[etaN]/.Jf[m_]:>2 Sinh[(qN+m)etaN]/(qN+m),25];
Ybar[l_,x_]:=Sqrt[(2 l+1)/(4 Pi)]LegendreP[l,x];
Kk[d_,l_,lp_,s_]:=2 Pi NIntegrate[Ybar[l,mu]Ybar[lp,(mu+s bb)/(1+s bb mu)]/(gN(1+s bb mu))^d,
   {mu,-1,1},WorkingPrecision->30,PrecisionGoal->18];
Dker[l_,lp_,ldd_]:=(1/gN)Kk[-1-qN,l,lp,-1]Kk[-qN,lp,ldd,1];

Print["T1  W-weight counts (l'=0..5), expect l'+1: ",Table[Length[W[lp]],{lp,0,5}]];
Print["T2  border tower from AXIS elements alone, D_{h,l',0} (h,l'<=5) vs direct: ",
  ToString[N@Max@Table[Abs[ev[Bord[h,lp]]-ev[Db[h,lp]]],{h,0,5},{lp,0,5}],InputForm],
  "   (0 => weights are outer-index independent)"];
Print["T3  MASTER formula vs kernel, ALL 216 elements: ",
  ToString[N@Max@Table[Abs[ev[Dmaster[l,lp,ldd]]-Dker[l,lp,ldd]],{l,0,5},{lp,0,5},{ldd,0,5}],InputForm]];

(* ---- (2) the normalised weights are ELEMENTARY (no J) ---- *)
cN[lp_]:=Table[red[W[lp][[j]] Jf[1-(j-1)]/(W[lp][[1]] Jf[1])],{j,1,Length[W[lp]]}];
Print["T4  normalised weights c^(l') are J-free (elementary in gamma,p): ",
  Table[FreeQ[cN[lp],Jf],{lp,1,5}]];

(* display them in gamma,p using gamma^2-p^2=1 (rationalise the gamma-denominator) *)
rg[e_]:=Expand[e]//.GAM^n_/;n>=2:>GAM^(n-2)(1+PP^2);
lauGP[L_]:=Module[{k=-Exponent[L,X,Min]},rg[Expand[Expand[L X^k]/.X->GAM+PP](GAM-PP)^k]];
gsplit[e_]:=Module[{ee=rg[e]},{ee/.GAM->0,Coefficient[ee,GAM]}];
canon[c_]:=Module[{t=Together[c],ne,no,de,do,rd},
  {ne,no}=gsplit[Numerator[t]]; {de,do}=gsplit[Denominator[t]];
  rd=Expand[de^2-(1+PP^2)do^2];
  Factor[Cancel[Expand[ne de-(1+PP^2)no do]/rd]]+GAM Factor[Cancel[Expand[no de-ne do]/rd]]];
fin[e_]:=Simplify[Cancel[Together[e/.(1+PP^2)->GAM^2]]];   (* gamma/(1+p^2) -> 1/gamma *)
cv[e_]:=fin[Module[{t=Together[e]},canon[Cancel[lauGP[Numerator[t]]/lauGP[Denominator[t]]]]]];
Print["T5  c^(1) = ",cv/@cN[1],"   (expect {1, -1/gamma})"];
Print["    c^(2)*(3+2p^2) = ",Expand[(3+2PP^2) cv[#]]&/@cN[2],
  "   (expect {3+2p^2, -6 gamma, 3} = published D020 bracket coefficients)"];
Print["    c^(3)*gamma(5+2p^2) = ",Expand[GAM (5+2PP^2) cv[#]]&/@cN[3]];
Print["T6  c^(2) coefficients match the published D020 bracket exactly: ",
  Simplify[Expand[(3+2PP^2) cv[#]]&/@cN[2] - {3+2PP^2,-6 GAM,3}]==={0,0,0}];

Print["--- general resummation harness done ---"];
