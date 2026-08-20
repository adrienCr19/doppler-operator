(* ::Package:: *)
(* 46_D555_substituted.wl — FULL SUBSTITUTION: unroll every relation to the bottom and put the
   result back into D_555, giving it in terms of the LOWEST index operator only, D_000.

   Unrolling (P0) then (P3) expresses each border element as a finite combination of the seed at
   integer-shifted spectral index,
        D_{h,l',0}(q) = Sum_m c^{(h,l')}_m(q) D_000(q+m),
   with c elementary (built from gamma, p and the brackets J_w).  Because the seed is
   reflection-symmetric, D_000(-3-q) = D_000(q) (N1 at l=l'=l''=0), the reflected border also
   collapses onto D_000 at shifted q:
        D_{h,l',0}(-3-q) = Sum_m c^{(h,l')}_m(-3-q) D_000(q-m).
   Hence, from (P2)+(N1),

     D_555(q) = [Sum_m c^{(5,5)}_m(q) D_000(q+m)] [Sum_m c^{(5,5)}_m(-3-q) D_000(q-m)]
                / [Sum_m c^{(0,5)}_m(q) D_000(q+m)] .

   Everything on the right is the single lowest operator D_000 at q, q+-1, ..., q+-5, times
   elementary coefficients.  The aberration kernel is used ONLY as the independent check. *)

g=(X+1/X)/2; p=(X-1/X)/2;                       (* X = e^eta *)
sh[e_,d_]:=e/.Jf[m_]:>Jf[m+d];
refl[e_]:=e/.Jf[m_]:>Jf[3-m];                   (* q -> -3-q *)
CC[l_]:=Sqrt[l^2/(4 l^2-1)]; red[e_]:=Cancel[Together[e]];
D000=Jf[2]Jf[1]/(4 g p^2);

(* ---------- 1. the ladders, symbolically ---------- *)
Da[0]=D000; Da[a_]/;a<0:=0;
Da[a_]:=Da[a]=red[(1/CC[a])((g/p)Da[a-1]-(1/p)(Jf[1]/Jf[2])sh[Da[a-1],1]-CC[a-1]Da[a-2])];
Tn[n_]:=(1/p^(n+1))Sum[Binomial[n,k](-g)^k Jf[n+2-k],{k,0,n}];
Jbr[a_]:=(-1)^a Sum[Coefficient[LegendreP[a,x],x,n]Tn[n],{n,0,a}];
aC[lp_]:=(g(2 lp+1)/(lp p)) sh[Jbr[lp],-1]/sh[Jbr[lp-1],-1];
bC[lp_]:=-((2 lp+1)/(lp p)) sh[Jbr[lp],-1]/sh[Jbr[lp-1],-2];
cC[lp_]:=If[lp>=2,-((lp-1)(2 lp+1)/(lp(2 lp-3))) sh[Jbr[lp],-1]/sh[Jbr[lp-2],-1],0];
Db[h_,0]:=Da[h]; Db[h_,lp_]/;lp<0:=0;
Db[h_,lp_]:=Db[h,lp]=red[aC[lp]Db[h,lp-1]+bC[lp]sh[Db[h,lp-1],-1]+cC[lp]Db[h,lp-2]];

(* ---------- 2. unroll onto D_000 shifts ---------- *)
(* axis:   Dhat_a(q)=Sum_k rho_k Dhat_0(q+k);  middle: D_{h,l',0}=J_l'(q-1) Sum_j U_j Dhat_h(q-j) *)
alp[l_]:=g(2 l+1)/(l p); bet[l_]:=-((2 l+1)/(l p)); kap[l_]:=-((l-1)(2 l+1)/(l(2 l-3)));
U[0]={p}; U[l_]/;l<0:={};
U[l_]:=U[l]=Module[{Um1=U[l-1],Um2=If[l>=2,U[l-2],{}],out},out=ConstantArray[0,l+1];
  Do[out[[j]]=out[[j]]+alp[l] Um1[[j]],{j,1,Length[Um1]}];
  Do[out[[j+1]]=out[[j+1]]+bet[l] Um1[[j]],{j,1,Length[Um1]}];
  If[l>=2,Do[out[[j]]=out[[j]]+kap[l] Um2[[j]],{j,1,Length[Um2]}]];Simplify/@out];
R[0]={1}; R[a_]/;a<0:={};
R[a_]:=R[a]=Module[{Rm1=R[a-1],Rm2=If[a>=2,R[a-2],{}],out},out=ConstantArray[0,a+1];
  Do[out[[k]]=out[[k]]+(1/CC[a])(g/p) Rm1[[k]],{k,1,Length[Rm1]}];
  Do[out[[k+1]]=out[[k+1]]-(1/CC[a])(1/p)(Jf[1]/Jf[2]) sh[Rm1[[k]],1],{k,1,Length[Rm1]}];
  If[a>=2,Do[out[[k]]=out[[k]]-(1/CC[a])CC[a-1] Rm2[[k]],{k,1,Length[Rm2]}]];red/@out];
rho[a_]:=rho[a]=Table[red[R[a][[k]] Jf[1+(k-1)]/Jf[1]],{k,1,Length[R[a]]}];
Dhat0[x_]:=sh[D000,x]/sh[Jf[1],x];
Dhat[a_,x_]:=Sum[rho[a][[k]] Dhat0[x+(k-1)],{k,1,Length[rho[a]]}];
border[h_,lp_]:=red[sh[Jbr[lp],-1] Sum[U[lp][[j]] Dhat[h,-(j-1)],{j,1,Length[U[lp]]}]];

(* extract the coefficient of D_000(q+m):  coef_m = J_l'(q-1) * A_m / J_{q+1+m} *)
cross[h_,lp_]:=Table[Sum[If[1<=j+m<=h+1,U[lp][[j]] rho[h][[j+m]],0],{j,1,lp+1}],{m,-lp,h}];
coef[h_,lp_]:=Table[red[sh[Jbr[lp],-1] cross[h,lp][[m+lp+1]]/sh[Jf[1],m]],{m,-lp,h}];
shiftsOf[h_,lp_]:=Range[-lp,h];

Print["T1  border as a combination of the SEED at shifted q:"];
Do[With[{h=c[[1]],lp=c[[2]]},
   Print["      D_{",h,",",lp,",0}(q) = Sum_m coef_m D_000(q+m),  m in ",shiftsOf[h,lp],
    "   exact: ",red[Sum[coef[h,lp][[i]] sh[D000,shiftsOf[h,lp][[i]]],{i,1,Length[coef[h,lp]]}]-Db[h,lp]]===0]],
  {c,{{1,0},{0,1},{1,1},{5,5},{0,5}}}];

(* ---------- 3. the seed is reflection-symmetric ---------- *)
Print["T2  seed reflection symmetry  D_000(-3-q) = D_000(q): ",red[refl[D000]-D000]===0];

(* ---------- 4. assemble D_555 purely from D_000 shifts ---------- *)
num1=Sum[coef[5,5][[i]] sh[D000,shiftsOf[5,5][[i]]],{i,1,Length[coef[5,5]]}];
num2=refl[num1];                                   (* q -> -3-q; seed shifts fold back by T2 *)
den =Sum[coef[0,5][[i]] sh[D000,shiftsOf[0,5][[i]]],{i,1,Length[coef[0,5]]}];
D555sub=red[num1 num2/den];
D555dir=red[Db[5,5] refl[Db[5,5]]/Db[0,5]];
Print["T3  D_555 assembled from D_000 shifts equals the direct product form: ",red[D555sub-D555dir]===0];

(* ---------- 5. independent numerical check against the kernel ---------- *)
bb=42/100; gN=1/Sqrt[1-bb^2]; etaN=ArcSinh[gN bb]; qN=37/100;
ev[e_]:=N[e/.X->Exp[etaN]/.Jf[m_]:>2 Sinh[(qN+m)etaN]/(qN+m),25];
Ybar[l_,x_]:=Sqrt[(2 l+1)/(4 Pi)]LegendreP[l,x];
Kk[d_,l_,lp_,s_]:=Kk[d,l,lp,s]=2 Pi NIntegrate[Ybar[l,mu]Ybar[lp,(mu+s bb)/(1+s bb mu)]/(gN(1+s bb mu))^d,
   {mu,-1,1},WorkingPrecision->30,PrecisionGoal->18];
Print["T4  vs aberration kernel: ",ToString[N[Abs[ev[D555sub]-(1/gN)Kk[-1-qN,5,5,-1]Kk[-qN,5,5,1]]],InputForm]];

(* ---------- 6. show the simplest fully-substituted cases explicitly ---------- *)
rg[e_]:=Expand[e]//.GAM^n_/;n>=2:>GAM^(n-2)(1+PP^2);
lauGP[L_]:=Module[{k=-Exponent[L,X,Min]},rg[Expand[Expand[L X^k]/.X->GAM+PP](GAM-PP)^k]];
cv[e_]:=Module[{t=Together[e]},Factor[Cancel[lauGP[Numerator[t]]/lauGP[Denominator[t]]]]];
Print["T5  simplest fully-substituted borders (coefficients of D_000(q+m)):"];
Do[With[{h=c[[1]],lp=c[[2]]},
   Print["      D_{",h,",",lp,",0}: shifts ",shiftsOf[h,lp]];
   Print["         coefficients ",cv/@coef[h,lp]]],
  {c,{{1,0},{0,1},{1,1}}}];
Print["T6  number of seed terms: D_{5,5,0} -> ",Length[coef[5,5]],
      " (shifts ",First[shiftsOf[5,5]],"..",Last[shiftsOf[5,5]],");  D_{0,5,0} -> ",Length[coef[0,5]]];
Print["--- D555 substitution harness done ---"];
