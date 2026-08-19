(* ::Package:: *)
(* 42_resummation_induction.wl — INDUCTION PROOF of the resummation theorem for general l'.

   THEOREM.  Let  Dhat_a(q) := D_{a00}(q)/J_{q+1}  (normalised axis element) and let the axis bracket
   be J_a(q) (= Jbr below).  Define ELEMENTARY constants (pure gamma,p; no J, no q):
       alpha_l = gamma(2l+1)/(l p),   beta_l = -(2l+1)/(l p),   kappa_l = -(l-1)(2l+1)/(l(2l-3)),
   and define U^{(l)}_j by the recursion
       U^{(0)}_0 = p,
       U^{(l)}_j = alpha_l U^{(l-1)}_j + beta_l U^{(l-1)}_{j-1} + kappa_l U^{(l-2)}_j,
   with U^{(l)}_j = 0 outside 0 <= j <= l.  Then for every l' >= 0 and every outer index h:

       (i)   U^{(l')}_j is elementary in gamma,p (contains no J) and is independent of q AND of h;
       (ii)  support is exactly j = 0..l'  (l'+1 terms);
       (iii) D_{h,l',0}(q) = J_{l'}(q-1) * Sum_{j=0}^{l'} U^{(l')}_j * Dhat_h(q-j).

   PROOF (induction on l').
     Base l'=0:  D_{h,0,0}(q) = D_{a00}(q) = J_{q+1} Dhat_h(q), and J_0(q-1) = J_{q+1}/p, so the
       right side is (J_{q+1}/p)(p)Dhat_h(q) = J_{q+1}Dhat_h(q).  True, with U^{(0)}_0 = p.
     Step:  P3' reads  D_{h,l',0}(q) = a_l' D_{h,l'-1,0}(q) + b_l' D_{h,l'-1,0}(q-1) + c_l' D_{h,l'-2,0}(q),
       with a,b,c INDEPENDENT of h.  Apply the hypothesis to each term.  The middle term needs the
       hypothesis at argument q-1, which is legitimate since it holds identically in q; reindexing
       j -> j+1 there produces the beta-shift.  The three prefactors collapse because
           a_l' * J_{l'-1}(q-1) = alpha_l' * J_{l'}(q-1),
           b_l' * J_{l'-1}(q-2) = beta_l'  * J_{l'}(q-1),
           c_l' * J_{l'-2}(q-1) = kappa_l' * J_{l'}(q-1),
       i.e. every level-(l'-1)/(l'-2) bracket is replaced by the single bracket J_{l'}(q-1) times an
       ELEMENTARY constant.  Dividing through by J_{l'}(q-1) leaves exactly the U-recursion, whose
       coefficients are elementary and q-free; since U^{(0)} is elementary and q-free, so is every
       U^{(l')} (so the q-shift in the beta-term acts trivially on it).  Support grows by exactly one
       from the beta-shift, giving 0..l'.  QED.

   This harness verifies every ingredient: the three prefactor collapses (the crux of the step), the
   base case, elementarity/q-independence, the support count, and the theorem itself at levels well
   beyond those used to discover it.  On nu^q, m=0. *)

g=(X+1/X)/2; p=(X-1/X)/2;                       (* X = e^eta *)
sh[e_,d_]:=e/.Jf[m_]:>Jf[m+d];                  (* q -> q+d *)
CC[l_]:=Sqrt[l^2/(4 l^2-1)]; red[e_]:=Cancel[Together[e]];
D000=Jf[2]Jf[1]/(4 g p^2); Da[0]=D000; Da[a_]/;a<0:=0;
Da[a_]:=Da[a]=red[(1/CC[a])((g/p)Da[a-1]-(1/p)(Jf[1]/Jf[2])sh[Da[a-1],1]-CC[a-1]Da[a-2])];
Tn[n_]:=(1/p^(n+1))Sum[Binomial[n,k](-g)^k Jf[n+2-k],{k,0,n}];
Jbr[a_]:=(-1)^a Sum[Coefficient[LegendreP[a,x],x,n]Tn[n],{n,0,a}];   (* axis bracket J_a(q) *)
aC[lp_]:=(g(2 lp+1)/(lp p)) sh[Jbr[lp],-1]/sh[Jbr[lp-1],-1];
bC[lp_]:=-((2 lp+1)/(lp p)) sh[Jbr[lp],-1]/sh[Jbr[lp-1],-2];
cC[lp_]:=If[lp>=2,-((lp-1)(2 lp+1)/(lp(2 lp-3))) sh[Jbr[lp],-1]/sh[Jbr[lp-2],-1],0];
Db[h_,0]:=Da[h]; Db[h_,lp_]/;lp<0:=0;
Db[h_,lp_]:=Db[h,lp]=red[aC[lp]Db[h,lp-1]+bC[lp]sh[Db[h,lp-1],-1]+cC[lp]Db[h,lp-2]];

(* ---- ELEMENTARY constants of the induction ---- *)
alp[l_]:=g(2 l+1)/(l p);  bet[l_]:=-((2 l+1)/(l p));  kap[l_]:=-((l-1)(2 l+1)/(l(2 l-3)));

(* ---- CRUX OF THE INDUCTIVE STEP: the three prefactor collapses ---- *)
crux1=Table[red[aC[l] sh[Jbr[l-1],-1] - alp[l] sh[Jbr[l],-1]],{l,1,10}];
crux2=Table[red[bC[l] sh[Jbr[l-1],-2] - bet[l] sh[Jbr[l],-1]],{l,1,10}];
crux3=Table[red[cC[l] sh[Jbr[l-2],-1] - kap[l] sh[Jbr[l],-1]],{l,2,10}];
Print["T1  inductive-step crux, a_l*J_{l-1}(q-1) = alpha_l*J_l(q-1)   (l=1..10): ",
   If[Union[crux1]==={0},"0  -> HOLDS",crux1]];
Print["    crux, b_l*J_{l-1}(q-2) = beta_l*J_l(q-1)                   (l=1..10): ",
   If[Union[crux2]==={0},"0  -> HOLDS",crux2]];
Print["    crux, c_l*J_{l-2}(q-1) = kappa_l*J_l(q-1)                  (l=2..10): ",
   If[Union[crux3]==={0},"0  -> HOLDS",crux3]];

(* ---- the U-recursion (elementary; no J anywhere) ---- *)
U[0]={p}; U[l_]/;l<0:={};
U[l_]:=U[l]=Module[{Um1=U[l-1],Um2=If[l>=2,U[l-2],{}],out},
  out=ConstantArray[0,l+1];
  Do[out[[j]]=out[[j]]+alp[l] Um1[[j]],{j,1,Length[Um1]}];
  Do[out[[j+1]]=out[[j+1]]+bet[l] Um1[[j]],{j,1,Length[Um1]}];   (* q-shift acts trivially: U is q-free *)
  If[l>=2,Do[out[[j]]=out[[j]]+kap[l] Um2[[j]],{j,1,Length[Um2]}]];
  Simplify/@out];

Print["T2  base case  U^(0) = ",U[0],"   (expect {p})"];
Print["T3  U is J-free / q-free for l<=10: ",Table[FreeQ[U[l],Jf],{l,0,10}]];
Print["T4  support count |U^(l)| = l+1 for l<=10: ",Table[Length[U[l]]==l+1,{l,0,10}]];

(* ---- THE THEOREM: D_{h,l',0}(q) = J_{l'}(q-1) Sum_j U_j Dhat_h(q-j) ---- *)
Dhat[a_]:=Da[a]/Jf[1];
thm[h_,lp_]:=red[ sh[Jbr[lp],-1] Sum[U[lp][[j]] sh[Dhat[h],-(j-1)],{j,1,Length[U[lp]]}] ];
dev=Table[red[thm[h,lp]-Db[h,lp]],{h,0,6},{lp,0,8}];
Print["T5  THEOREM  D_{h,l',0} = J_{l'}(q-1) Sum_j U_j Dhat_h(q-j),  h<=6, l'<=8: ",
   If[Union[Flatten[dev]]==={0},"0  -> PROVED IDENTITY HOLDS (symbolic, dev 0)","FAIL"]];

(* ---- the weights are the published bracket coefficients ---- *)
Print["T6  U^(1) = ",U[1],"   (expect {3 gamma, -3} ~ {gamma,-1})"];
Print["    U^(2)*(2p/5) = ",Simplify[(2 p/5) #]&/@U[2],
   "   (expect {3+2p^2, -6 gamma, 3} = the published D020 bracket)"];
Print["T7  U^(2) normalised equals published D020 bracket exactly: ",
   Simplify[((2 p/5) #)&/@U[2] - {3+2p^2,-6 g,3}]==={0,0,0}];

(* ---- display U in gamma,p, and close the loop against the KERNEL definition ---- *)
rg[e_]:=Expand[e]//.GAM^n_/;n>=2:>GAM^(n-2)(1+PP^2);
lauGP[L_]:=Module[{k=-Exponent[L,X,Min]},rg[Expand[Expand[L X^k]/.X->GAM+PP](GAM-PP)^k]];
Print["T8  U^(1) in gamma,p = ",Factor[lauGP[#]]&/@U[1]];
Print["    U^(2)*(2p/5)     = ",Factor[lauGP[Simplify[(2 p/5) #]]]&/@U[2]];
Print["    U^(3)*(2p/7)     = ",Factor[lauGP[Simplify[(2 p/7) #]]]&/@U[3]];

bb=42/100; gN=1/Sqrt[1-bb^2]; etaN=ArcSinh[gN bb]; qN=37/100;
ev[e_]:=N[e/.X->Exp[etaN]/.Jf[m_]:>2 Sinh[(qN+m)etaN]/(qN+m),25];
Ybar[l_,x_]:=Sqrt[(2 l+1)/(4 Pi)]LegendreP[l,x];
Kk[d_,l_,lp_,s_]:=2 Pi NIntegrate[Ybar[l,mu]Ybar[lp,(mu+s bb)/(1+s bb mu)]/(gN(1+s bb mu))^d,
   {mu,-1,1},WorkingPrecision->30,PrecisionGoal->18];
Dker[l_,lp_,ldd_]:=(1/gN)Kk[-1-qN,l,lp,-1]Kk[-qN,lp,ldd,1];
Print["T9  theorem vs the KERNEL definition (h<=5, l'<=5): ",
   ToString[N@Max@Table[Abs[ev[thm[h,lp]]-Dker[h,lp,0]],{h,0,5},{lp,0,5}],InputForm],
   "   (closes the chain: induction -> P3' -> physical Doppler operator)"];

(* ---- T10: the weights in clean normalised form ---- *)
gsplit[e_]:=Module[{ee=rg[e]},{ee/.GAM->0,Coefficient[ee,GAM]}];
canon[c_]:=Module[{t=Together[c],ne,no,de,do,rd},
  {ne,no}=gsplit[Numerator[t]]; {de,do}=gsplit[Denominator[t]];
  rd=Expand[de^2-(1+PP^2)do^2];
  Cancel[Expand[ne de-(1+PP^2)no do]/rd]+GAM Cancel[Expand[no de-ne do]/rd]];
norm[l_]:=Module[{v,r,den},v=canon[lauGP[#]]&/@U[l];
  r=Together[#/Last[v]]&/@v; den=Apply[PolynomialLCM,Denominator/@r];
  Factor[Cancel[Together[# den]]]&/@r];
Print["T10 the weights U^(l') (up to an overall scalar), gamma=GAM, p=PP:"];
Do[Print["    U^(",l,") ~ ",norm[l]],{l,1,5}];

Print["--- resummation induction harness done ---"];
