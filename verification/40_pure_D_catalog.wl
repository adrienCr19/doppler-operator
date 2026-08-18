(* ::Package:: *)
(* 40_pure_D_catalog.wl — Generate CLOSED FORMS for every Doppler operator D_{l l' l''}(q) up to
   D_{555} using ONLY the pure-D recurrence ladder (P0 axis, P3' middle, N1 reflection, P2 fill).
   The aberration kernel is NEVER used to generate anything; it appears only as an independent
   numeric check.  Exports data/closed_forms_pureD.json.  On nu^q, m=0. *)
g=(X+1/X)/2; p=(X-1/X)/2;                       (* X = e^eta, internal only *)
sh[e_,d_]:=e/.Jf[m_]:>Jf[m+d];                  (* q -> q+d *)
refl[e_]:=e/.Jf[m_]:>Jf[3-m];                   (* N1: q -> -3-q  (J-index m -> 3-m) *)
CC[l_]:=Sqrt[l^2/(4 l^2-1)]; red[e_]:=Cancel[Together[e]];
(* seed *)
D000=Jf[2]Jf[1]/(4 g p^2);
(* P0: outer axis ladder *)
Da[0]=D000; Da[a_]/;a<0:=0;
Da[a_]:=Da[a]=red[(1/CC[a])((g/p)Da[a-1]-(1/p)(Jf[1]/Jf[2])sh[Da[a-1],1]-CC[a-1]Da[a-2])];
(* P3': middle ladder, elementary J-bracket coefficients *)
Tn[n_]:=(1/p^(n+1))Sum[Binomial[n,k](-g)^k Jf[n+2-k],{k,0,n}];
Jbr[a_]:=(-1)^a Sum[Coefficient[LegendreP[a,x],x,n]Tn[n],{n,0,a}];
aC[lp_]:=(g(2 lp+1)/(lp p)) sh[Jbr[lp],-1]/sh[Jbr[lp-1],-1];
bC[lp_]:=-((2 lp+1)/(lp p)) sh[Jbr[lp],-1]/sh[Jbr[lp-1],-2];
cC[lp_]:=If[lp>=2,-((lp-1)(2 lp+1)/(lp(2 lp-3))) sh[Jbr[lp],-1]/sh[Jbr[lp-2],-1],0];
Db[h_,0]:=Da[h]; Db[h_,lp_]/;lp<0:=0;
Db[h_,lp_]:=Db[h,lp]=red[aC[lp]Db[h,lp-1]+bC[lp]sh[Db[h,lp-1],-1]+cC[lp]Db[h,lp-2]];
(* P2 + N1: any element from the border tower *)
Dfull[l_,lp_,ldd_]:=red[Db[l,lp] refl[Db[ldd,lp]]/Db[0,lp]];
(* ---- display converter: X-scalar -> polynomial in GAM(gamma),PP(p); reduce GAM^2->1+PP^2 ---- *)
(* split polynomial in GAM (using gamma^2=1+p^2) into {A,B} with expr = A + GAM B, A,B in PP *)
gsplit[e_]:=Module[{ee=Expand[e]//.GAM^n_/;n>=2:>(1+PP^2)^Quotient[n,2] GAM^Mod[n,2]},
   {ee/.GAM->0, Coefficient[ee,GAM]}];
(* canonical rationalized form  A(p) + gamma B(p) *)
canon[c_]:=Module[{t=Together[c],ne,no,de,do,rd},
   {ne,no}=gsplit[Numerator[t]]; {de,do}=gsplit[Denominator[t]];
   rd=Expand[de^2-(1+PP^2)do^2];
   Cancel[Expand[ne de-(1+PP^2)no do]/rd]+GAM Cancel[Expand[no de-ne do]/rd]];
sgp[c_]:=Module[{r=ExpToTrig[c/.X->Exp[et]]},
  r=r//.{Csch[x_]:>1/Sinh[x],Sech[x_]:>1/Cosh[x],Coth[x_]:>Cosh[x]/Sinh[x],Tanh[x_]:>Sinh[x]/Cosh[x]};
  r=r/.{Cosh[k_. et]:>ChebyshevT[k,GAM],Sinh[k_. et]:>PP ChebyshevU[k-1,GAM]};
  (* canonical A(p)+gamma B(p), then restore site/paper style: (1+p^2)->gamma^2, cancel *)
  Cancel[Together[Factor[canon[r]]/.(1+PP^2)->GAM^2]]];
vars=Table[Jf[i],{i,-4,10}];
toJ[e_]:=Module[{t=Together[e],num,den,cr},den=Denominator[t];num=Numerator[t];cr=CoefficientRules[num,vars];
  Total[(sgp[Cancel[#[[2]]/den]] Times@@(vars^#[[1]]))&/@cr]];
(* serialize to LaTeX with GAM->gamma, PP->p at the STRING level (no symbol leak) *)
texit[Dc_]:=StringReplace[ToString[TeXForm[Dc/.Jf[m_]:>Subscript[JJJ,m]]],
   {"\\text{GAM}"->"\\gamma","GAM"->"\\gamma","\\text{PP}"->"p","PP"->"p","\\text{JJJ}"->"J","JJJ"->"J"}];
(* ---- kernel reference (CHECK ONLY) ---- *)
bb=42/100; gN=1/Sqrt[1-bb^2]; etaN=ArcSinh[gN bb]; qN=37/100;
Jq[m_]:=2 Sinh[(qN+m)etaN]/(qN+m);
evX[e_]:=N[e/.X->Exp[etaN]/.Jf[m_]:>Jq[m],25];                (* eval X-form *)
evD[e_]:=N[e/.{GAM->gN,PP->gN bb}/.Jf[m_]:>Jq[m],25];        (* eval converted form *)
Ybar[l_,x_]:=Sqrt[(2 l+1)/(4 Pi)]LegendreP[l,x];
Kk[d_,l_,lp_,s_]:=2 Pi NIntegrate[Ybar[l,mu]Ybar[lp,(mu+s bb)/(1+s bb mu)]/(gN(1+s bb mu))^d,{mu,-1,1},WorkingPrecision->30,PrecisionGoal->18];
Dker[l_,lp_,ldd_]:=(1/gN)Kk[-1-qN,l,lp,-1]Kk[-qN,lp,ldd,1];
(* ---- run: verify ladder + converter for all 216, export ---- *)
maxlad=0; maxconv=0; recs={};
Do[Module[{Dx,Dc,vk},Dx=Dfull[l,lp,ldd]; vk=Dker[l,lp,ldd];
   maxlad=Max[maxlad,Abs[evX[Dx]-vk]];
   Dc=toJ[Dx]; maxconv=Max[maxconv,Abs[evD[Dc]-vk]];
   recs=Append[recs,<|"l"->l,"lp"->lp,"ldd"->ldd,"form"->texit[Dc]|>]],
  {l,0,5},{lp,0,5},{ldd,0,5}];
Print["ladder vs kernel  (216 elts): ",N[maxlad]];
Print["converted-form vs kernel:     ",N[maxconv]];
Export["/Users/adriencristian/Desktop/rSZ/doppler-operator/data/closed_forms_pureD.json",recs];
Print["exported ",Length[recs]," closed forms to data/closed_forms_pureD.json"];
(* print a few clean ones for the log *)
Print["D000: ",texit[toJ[Da[0]]]];
Print["D100: ",texit[toJ[Da[1]]]];
Print["D020: ",texit[toJ[Db[0,2]]]];
Print["D110: ",texit[toJ[Db[1,1]]]];
Print["--- pure-D catalog harness done ---"];
