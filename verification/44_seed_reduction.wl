(* ::Package:: *)
(* 44_seed_reduction.wl — ALL THE WAY DOWN: every Doppler operator reduced to the SEED alone.

   Both ladders resum in the normalised variable  Dhat_a(q) := D_{a00}(q)/J_{q+1} :

     AXIS   (P0 unrolled):   Dhat_a(q) = Sum_{k=0}^{a} rho^{(a)}_k Dhat_0(q+k)
     MIDDLE (P3' unrolled):  D_{h,l',0}(q) = J_{l'}(q-1) Sum_{j=0}^{l'} U^{(l')}_j Dhat_h(q-j)

   with BOTH weight families rho^{(a)}_k and U^{(l')}_j ELEMENTARY in gamma,p (no J, no q).
   Composing them, and closing with the rank-one fill (P2) + reflection (N1), every element is an
   explicit finite combination of the SEED  Dhat_0(x) = D_{000}(x)/J_{x+1}  at integer-shifted x:

     D_{h,l',0}(q) = J_{l'}(q-1) Sum_{j,k} U^{(l')}_j rho^{(h)}_k Dhat_0(q-j+k)
     D_{l l' l''}(q) = D_{l l' 0}(q) D_{l'' l' 0}(-3-q) / D_{0 l' 0}(q)

   and the common prefactor J_{l'}(q-1) CANCELS in the last line, leaving for the top element

     D_555(q) = J_5(-4-q) * [Sum_{j,k} U_j rho_k Dhat_0(q-j+k)]
                          * [Sum_{j,k} U_j rho_k Dhat_0(-3-q-j+k)]
                          / [Sum_m U_m Dhat_0(q-m)] .

   On nu^q, m=0.  (Operator form: q -> -Ohat.) *)

g=(X+1/X)/2; p=(X-1/X)/2;
sh[e_,d_]:=e/.Jf[m_]:>Jf[m+d];
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
Dfull[l_,lp_,ldd_]:=red[Db[l,lp] (Db[ldd,lp]/.Jf[m_]:>Jf[3-m])/Db[0,lp]];

(* ---- middle weights U (elementary; harness 42) ---- *)
alp[l_]:=g(2 l+1)/(l p); bet[l_]:=-((2 l+1)/(l p)); kap[l_]:=-((l-1)(2 l+1)/(l(2 l-3)));
U[0]={p}; U[l_]/;l<0:={};
U[l_]:=U[l]=Module[{Um1=U[l-1],Um2=If[l>=2,U[l-2],{}],out},out=ConstantArray[0,l+1];
  Do[out[[j]]=out[[j]]+alp[l] Um1[[j]],{j,1,Length[Um1]}];
  Do[out[[j+1]]=out[[j+1]]+bet[l] Um1[[j]],{j,1,Length[Um1]}];
  If[l>=2,Do[out[[j]]=out[[j]]+kap[l] Um2[[j]],{j,1,Length[Um2]}]];Simplify/@out];

(* ---- AXIS: unroll P0 into shift coefficients R, then normalise to rho ---- *)
R[0]={1}; R[a_]/;a<0:={};
R[a_]:=R[a]=Module[{Rm1=R[a-1],Rm2=If[a>=2,R[a-2],{}],out},out=ConstantArray[0,a+1];
  Do[out[[k]]=out[[k]]+(1/CC[a])(g/p) Rm1[[k]],{k,1,Length[Rm1]}];
  Do[out[[k+1]]=out[[k+1]]-(1/CC[a])(1/p)(Jf[1]/Jf[2]) sh[Rm1[[k]],1],{k,1,Length[Rm1]}];
  If[a>=2,Do[out[[k]]=out[[k]]-(1/CC[a])CC[a-1] Rm2[[k]],{k,1,Length[Rm2]}]];
  red/@out];
rho[a_]:=rho[a]=Table[red[R[a][[k]] Jf[1+(k-1)]/Jf[1]],{k,1,Length[R[a]]}];

Print["T1  axis unrolled: D_{a00}(q) = Sum_k R^(a)_k D_000(q+k), a=0..5: ",
  If[Union[Table[red[Sum[R[a][[k]] sh[D000,(k-1)],{k,1,Length[R[a]]}]-Da[a]],{a,0,5}]]==={0},
     "0 -> EXACT","MISMATCH"]];
Print["T2  normalised axis weights rho^(a) are J-free (elementary in gamma,p), a=0..8: ",
  Table[FreeQ[rho[a],Jf],{a,0,8}]];
Print["T3  support |rho^(a)| = a+1, a=0..8: ",Table[Length[rho[a]]==a+1,{a,0,8}]];

(* ---- seed-only reduction ---- *)
Dhat0[x_]:=sh[D000,x]/sh[Jf[1],x];                       (* Dhat_0(q+x) *)
Dhat[a_,x_]:=Sum[rho[a][[k]] Dhat0[x+(k-1)],{k,1,Length[rho[a]]}];
border[h_,lp_]:=red[ sh[Jbr[lp],-1] Sum[U[lp][[j]] Dhat[h,-(j-1)],{j,1,Length[U[lp]]}] ];
Print["T4  border tower from the SEED alone, D_{h,l',0} (h<=5, l'<=5): ",
  If[Union[Flatten[Table[red[border[h,lp]-Db[h,lp]],{h,0,5},{lp,0,5}]]]==={0},
     "0 -> EXACT (symbolic)","MISMATCH"]];

full[l_,lp_,ldd_]:=red[ border[l,lp] (border[ldd,lp]/.Jf[m_]:>Jf[3-m]) / border[0,lp] ];
Print["T5  FULL element from the seed alone, all l,l',l''<=3: ",
  If[Union[Flatten[Table[red[full[l,lp,ldd]-Dfull[l,lp,ldd]],{l,0,3},{lp,0,3},{ldd,0,3}]]]==={0},
     "0 -> EXACT (symbolic)","MISMATCH"]];
Print["T6  D_555 from the seed alone: ",
  If[red[full[5,5,5]-Dfull[5,5,5]]===0,"0 -> EXACT (symbolic)","MISMATCH"]];

(* ---- display the two elementary weight families ---- *)
rg[e_]:=Expand[e]//.GAM^n_/;n>=2:>GAM^(n-2)(1+PP^2);
lauGP[L_]:=Module[{k=-Exponent[L,X,Min]},rg[Expand[Expand[L X^k]/.X->GAM+PP](GAM-PP)^k]];
gsplit[e_]:=Module[{ee=rg[e]},{ee/.GAM->0,Coefficient[ee,GAM]}];
canon[c_]:=Module[{t=Together[c],ne,no,de,do,rd},{ne,no}=gsplit[Numerator[t]];{de,do}=gsplit[Denominator[t]];
  rd=Expand[de^2-(1+PP^2)do^2];
  Factor[Cancel[Expand[ne de-(1+PP^2)no do]/rd]]+GAM Factor[Cancel[Expand[no de-ne do]/rd]]];
nrm[v_]:=Module[{c=canon[lauGP[#]]&/@v,r,den},r=Together[#/Last[c]]&/@c;
  den=Apply[PolynomialLCM,Denominator/@r]; Factor[Cancel[Together[# den]]]&/@r];
Print["T7  axis weights rho^(a) (up to an overall scalar):"];
Do[Print["      rho^(",a,") ~ ",nrm[rho[a]]],{a,1,4}];
Print["    middle weights U^(l') (up to an overall scalar):"];
Do[Print["      U^(",l,") ~ ",nrm[U[l]]],{l,1,4}];
(* ---- T8: the two weight families are PROPORTIONAL; the double sum is a cross-correlation ---- *)
prop[n_]:=Module[{r=Simplify[rho[n]/U[n]]},{Length[Union[r]]==1,Simplify[First[r]]}];
Print["T8  rho^(n) proportional to U^(n)?  {proportional, ratio}, n=1..6:"];
Do[Print["      n=",n,": ",prop[n]],{n,1,6}];
(* border = J_l'(q-1) Sum_{j,k} U^{(l')}_j rho^{(h)}_k Dhat_0(q+k-j) = J Sum_m A_m Dhat_0(q+m),
   A_m = Sum_j U^{(l')}_j rho^{(h)}_{j+m}   (cross-correlation of the two weight vectors)      *)
cross[h_,lp_]:=Table[Sum[If[1<=j+m<=h+1,U[lp][[j]] rho[h][[j+m]],0],{j,1,lp+1}],{m,-lp,h}];
borderX[h_,lp_]:=red[sh[Jbr[lp],-1] Sum[cross[h,lp][[m+lp+1]] Dhat0[m],{m,-lp,h}]];
Print["T9  cross-correlation form  D_{h,l',0}(q) = J_{l'}(q-1) Sum_m A_m Dhat_0(q+m), h,l'<=5: ",
  If[Union[Flatten[Table[red[borderX[h,lp]-Db[h,lp]],{h,0,5},{lp,0,5}]]]==={0},
     "0 -> EXACT (symbolic)","MISMATCH"]];
Print["    for h=l'=5 the shift range is m = -5..5 (11 terms) vs 36 in the double sum"];

Print["--- seed-reduction harness done ---"];
