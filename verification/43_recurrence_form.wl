(* ::Package:: *)
(* 43_recurrence_form.wl — Every Doppler operator written in terms of OTHER Doppler operators.
   For each (l,l',l'') exactly one pure-D one-step rule applies:

     (P2+N1)  l''>0 :  D_{l l' l''}(q) = D_{l l' 0}(q) D_{l'' l' 0}(-3-q) / D_{0 l' 0}(q)
     (P3)     l''=0, l'>=1 :
        D_{h l' 0}(q) = (1/C_l') [ (D_{l',00}(q-1)/p)( g D_{h,l'-1,0}(q)/D_{l'-1,00}(q-1)
                                     - (J_{q-1}/J_q) D_{h,l'-1,0}(q-1)/D_{l'-1,00}(q-2) )
                                  - C_{l'-1} (D_{l',00}(q-1)/D_{l'-2,00}(q-1)) D_{h,l'-2,0}(q) ]
     (P0)     l''=0, l'=0, l>=1 :
        D_{a00}(q) = (1/C_a)[ (g/p) D_{a-1,00}(q) - (1/p)(J_{1+q}/J_{2+q}) D_{a-1,00}(q+1)
                              - C_{a-1} D_{a-2,00}(q) ]
     (seed)   (0,0,0) : D_000(q) = J_{q+1}J_{q+2}/(4 g p^2)

   Every operator on the right is a physical Doppler operator (at a shifted spectral index);
   the only scalars are gamma, p and the elementary bracket ratios.  This harness checks each
   rule, element by element, against the aberration-kernel definition (used ONLY as the check),
   and exports the applicable rule + parents for every element.  On nu^q, m=0. *)

prec=30;
Ybar[l_,x_]:=Sqrt[(2 l+1)/(4 Pi)]LegendreP[l,x];
bb=42/100; gN=1/Sqrt[1-bb^2]; pN=gN bb; etaN=ArcSinh[pN];
Kk[d_?NumericQ,l_,lp_,s_]:=Kk[d,l,lp,s]=2 Pi NIntegrate[
   Ybar[l,mu]Ybar[lp,(mu+s bb)/(1+s bb mu)]/(gN(1+s bb mu))^d,{mu,-1,1},
   WorkingPrecision->prec,PrecisionGoal->18];
Kk[d_,l_,lp_,s_]/;l<0||lp<0:=0;
Dop[l_,lp_,ldd_,q_]/;l<0||lp<0||ldd<0:=0;
Dop[l_,lp_,ldd_,q_]:=(1/gN)Kk[-1-q,l,lp,-1]Kk[-q,lp,ldd,1];      (* reference only *)
JJ[w_]:=2 Sinh[w etaN]/w;
CC[l_]:=If[l<=0,0,Sqrt[l^2/(4 l^2-1)]];
rep[name_,dev_]:=Print[name,": ",If[TrueQ[N[dev]<10^-15],"PASS","FAIL"],"  (worst dev = ",ToString[N[dev],InputForm],")"];

(* ---- (P0) axis ladder: D_{a00} from D_{a-1,00}, D_{a-2,00} ---- *)
p0[a_,q_]:=(1/CC[a])((gN/pN)Dop[a-1,0,0,q]-(1/pN)(JJ[1+q]/JJ[2+q])Dop[a-1,0,0,q+1]-CC[a-1]Dop[a-2,0,0,q]);
rep["P0  axis, D_{a00} from one index down (a=1..5)",
  Max@Table[Abs[p0[a,q]-Dop[a,0,0,q]],{a,1,5},{q,{37/100,-(13/10)}}]];

(* ---- (P3) middle step: D_{h,l',0} from D_{h,l'-1,0}, D_{h,l'-2,0} ---- *)
p3[h_,lp_,q_]:=(1/CC[lp])(
   (Dop[lp,0,0,q-1]/pN)(gN Dop[h,lp-1,0,q]/Dop[lp-1,0,0,q-1]
      -(JJ[q-1]/JJ[q]) Dop[h,lp-1,0,q-1]/Dop[lp-1,0,0,q-2])
   -If[lp>=2,CC[lp-1](Dop[lp,0,0,q-1]/Dop[lp-2,0,0,q-1])Dop[h,lp-2,0,q],0]);
rep["P3  middle, D_{h,l',0} from one middle-index down (h<=5, l'=1..5)",
  Max@Table[Abs[p3[h,lp,q]-Dop[h,lp,0,q]],{h,0,5},{lp,1,5},{q,{37/100,-(13/10)}}]];

(* ---- (P2+N1) interior fill: D_{l l' l''} from the two borders ---- *)
p2[l_,lp_,ldd_,q_]:=Dop[l,lp,0,q] Dop[ldd,lp,0,-3-q]/Dop[0,lp,0,q];
rep["P2+N1  interior, D_{l l' l''} from border operators (all l,l'<=5, l''=1..5)",
  Max@Table[Abs[p2[l,lp,ldd,q]-Dop[l,lp,ldd,q]],{l,0,5},{lp,0,5},{ldd,1,5},{q,{37/100}}]];

(* ---- coverage: exactly one rule applies to each of the 216 elements ---- *)
rule[l_,lp_,ldd_]:=Which[ldd>0,"P2+N1", lp>=1,"P3", l>=1,"P0", True,"seed"];
cnt=Counts[Flatten[Table[rule[l,lp,ldd],{l,0,5},{lp,0,5},{ldd,0,5}]]];
Print["coverage of the 216 elements by rule: ",cnt,"   total = ",Total[Values[cnt]]];

(* ---- export the reduction table (rule + parent operators) ---- *)
parents[l_,lp_,ldd_]:=Which[
  ldd>0, {{l,lp,0,"q"},{ldd,lp,0,"-3-q"},{0,lp,0,"q"}},
  lp>=1, Join[{{lp,0,0,"q-1"},{lp-1,0,0,"q-1"},{lp-1,0,0,"q-2"},{l,lp-1,0,"q"},{l,lp-1,0,"q-1"}},
              If[lp>=2,{{lp-2,0,0,"q-1"},{l,lp-2,0,"q"}},{}]],
  l>=1, {{l-1,0,0,"q"},{l-1,0,0,"q+1"},{l-2,0,0,"q"}},
  True, {}];
recs=Flatten[Table[<|"l"->l,"lp"->lp,"ldd"->ldd,"rule"->rule[l,lp,ldd],
   "parents"->parents[l,lp,ldd]|>,{l,0,5},{lp,0,5},{ldd,0,5}],2];
Export["/Users/adriencristian/Desktop/rSZ/doppler-operator/data/recurrence_form.json",recs];
Print["exported ",Length[recs]," one-step reductions to data/recurrence_form.json"];
Print["example — D555: ",rule[5,5,5]," with parents ",parents[5,5,5]];
Print["--- recurrence-form harness done ---"];
