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
(* ---- FACTORED closed form:  D = S(gamma,p) * B1(J) * B2(J), two linear-in-J brackets ----
   Every element is bilinear in the J's (rank-one structure of P2), so it factors into exactly two
   brackets.  The factorization only exists modulo gamma^2 = 1+p^2, so we factor in the rational
   X = e^eta representation (where that relation is automatic) and convert back.                    *)
JV=Table[Jf[i],{i,-8,14}];  hasJ[e_]:=!FreeQ[e,Jf];
rg[e_]:=Expand[e]//.GAM^n_/;n>=2:>GAM^(n-2)(1+PP^2);          (* gamma^2 -> 1+p^2 *)
(* Laurent polynomial in X -> (gamma,p), denominator-free:  X -> GAM+PP,  1/X -> GAM-PP *)
lauGP[L_]:=Module[{k=-Exponent[L,X,Min]},rg[Expand[Expand[L X^k]/.X->GAM+PP](GAM-PP)^k]];
numContent[e_]:=Module[{cr=CoefficientRules[Expand[e],Join[JV,{GAM,PP}]]},
   If[cr==={},1,Apply[GCD,cr[[All,2]]]]];
simpS[s_]:=Module[{r=Simplify[Together[s]/.GAM->Sqrt[1+PP^2],PP>0]},
   Simplify[r/.Sqrt[1+PP^2]->GAM/.(1+PP^2)->GAM^2,GAM>0&&PP>0]];
facD[l_,lp_,ldd_]:=Module[{D0,fl,br,cen,Bc,Bgp,ks,Bn,Sx},
  D0=Together[Dfull[l,lp,ldd]];
  fl=FactorList[D0];
  br=Flatten[ConstantArray[#[[1]],#[[2]]]&/@Select[fl,hasJ[#[[1]]]&]];
  cen=(Exponent[#,X]+Exponent[#,X,Min])/2&/@br;               (* centre each bracket in X *)
  Bc=MapThread[Expand[#1 X^(-#2)]&,{br,cen}];
  Bgp=lauGP/@Bc; ks=numContent/@Bgp; Bn=MapThread[Expand[#1/#2]&,{Bgp,ks}];
  Sx=Cancel[Together[D0/(Times@@br)]];                         (* J-free scalar *)
  {simpS[lauGP[Sx] Times@@ks (GAM+PP)^Total[cen]], Collect[#,JV,Factor]&/@Bn}];
(* serialize S * B1 * B2 to LaTeX as  S \left(B1\right)\left(B2\right) *)
tx[e_]:=StringReplace[ToString[TeXForm[e/.Jf[m_]:>Subscript[JJJ,m]]],
  {"\\text{GAM}"->"\\gamma","GAM"->"\\gamma","\\text{PP}"->"p","PP"->"p","\\text{JJJ}"->"J","JJJ"->"J"}];
texit[f_]:=StringJoin[tx[f[[1]]],StringJoin[("\\left("<>tx[#]<>"\\right)")&/@f[[2]]]];

(* ---- kernel reference (CHECK ONLY) ---- *)
bb=42/100; gN=1/Sqrt[1-bb^2]; etaN=ArcSinh[gN bb]; qN=37/100;
Jq[m_]:=2 Sinh[(qN+m)etaN]/(qN+m);
evX[e_]:=N[e/.X->Exp[etaN]/.Jf[m_]:>Jq[m],25];                (* eval X-form *)
evD[e_]:=N[e/.{GAM->gN,PP->gN bb}/.Jf[m_]:>Jq[m],25];        (* eval converted form *)
Ybar[l_,x_]:=Sqrt[(2 l+1)/(4 Pi)]LegendreP[l,x];
Kk[d_,l_,lp_,s_]:=2 Pi NIntegrate[Ybar[l,mu]Ybar[lp,(mu+s bb)/(1+s bb mu)]/(gN(1+s bb mu))^d,{mu,-1,1},WorkingPrecision->30,PrecisionGoal->18];
Dker[l_,lp_,ldd_]:=(1/gN)Kk[-1-qN,l,lp,-1]Kk[-qN,lp,ldd,1];
(* ---- run: verify ladder + factored forms for all 216, export ---- *)
maxlad=0; maxfac=0; recs={};
Do[Module[{Dx,f,vk},Dx=Dfull[l,lp,ldd]; vk=Dker[l,lp,ldd];
   maxlad=Max[maxlad,Abs[evX[Dx]-vk]];                       (* pure-D ladder vs kernel *)
   f=facD[l,lp,ldd];
   maxfac=Max[maxfac,Abs[evD[f[[1]] Times@@f[[2]]]-vk]];     (* factored form vs kernel *)
   recs=Append[recs,<|"l"->l,"lp"->lp,"ldd"->ldd,"form"->texit[f],"nbrackets"->Length[f[[2]]]|>]],
  {l,0,5},{lp,0,5},{ldd,0,5}];
Print["ladder vs kernel     (216 elts): ",N[maxlad]];
Print["FACTORED form vs kernel:         ",N[maxfac]];
Print["all elements factor into 2 brackets: ",
   If[Union[#["nbrackets"]&/@recs]==={2}||Union[#["nbrackets"]&/@recs]==={1,2},
      Union[#["nbrackets"]&/@recs],"UNEXPECTED"]];
Export["/Users/adriencristian/Desktop/rSZ/doppler-operator/data/closed_forms_pureD.json",recs];
Print["exported ",Length[recs]," factored closed forms to data/closed_forms_pureD.json"];
(* a few, to compare against the published expressions *)
Print["D000: ",texit[facD[0,0,0]]];
Print["D101: ",texit[facD[1,0,1]],"   [published: 3/(4 g p^4)(g J_{1+q}-J_q)(g J_{2+q}-J_{3+q})]"];
Print["D020: ",texit[facD[0,2,0]],"   [published N8 Eq.1.3]"];
Print["D121: ",texit[facD[1,2,1]],"   [published N16 box: 15/(16 g p^8) W(q) W(-3-q)]"];
Print["--- pure-D catalog harness done ---"];
