(* ::Package:: *)
(* 26_all_index_ode.wl — a rapidity relation of D moving ALL THREE indices.
   N3 = d/deta as OUTER shifts (l,ldd);  mirror G5b = same d/deta as MIDDLE shifts (l1,l2).
   Both equal d/deta X, so: (i) they are EQUAL -> the algebraic CROSS-SLOT BALANCE (all three indices, no deriv);
   (ii) 1/2 N3 + 1/2 G5b = d/deta X moves all three at once.  Physical D = diagonal l1=l2=l'. m=0, q=1/2. *)
prec=40; nintg[f_]:=NIntegrate[f,{mu,-1,1},WorkingPrecision->prec,PrecisionGoal->25,MaxRecursion->25];
Ybar[l_,x_]:=Sqrt[(2l+1)/(4Pi)]LegendreP[l,x]; ga[b_]:=1/Sqrt[1-b^2];
KmR[d_?NumericQ,l_,lp_,b_]:=Module[{gg=ga[b]},2Pi nintg[Ybar[l,mu]Ybar[lp,(mu+b)/(1+b mu)]/(gg(1+b mu))^d]];
Km[d_?NumericQ,l_,lp_,b_]/;l<0||lp<0:=0; Km[d_?NumericQ,l_,lp_,b_]:=Km[d,l,lp,b]=KmR[d,l,lp,b];
Kp[d_?NumericQ,l_,lp_,b_]:=Km[d,l,lp,-b]; Cm0[l_]:=If[l<=0,0,Sqrt[l^2/(4l^2-1)]];
report[name_,dev_,tol_:10^-6]:=Print[name,": ",If[TrueQ[N[dev]<tol],"PASS","FAIL"],"  (dev = ",ScientificForm[N[dev],3],")"];
qv=1/2; bb=3/10; eta0=ArcTanh[bb]; de=1/1000000;
X[l_,l1_,l2_,ldd_,b_]:=(1/ga[b])Kp[-1-qv,l,l1,b]Km[-qv,l2,ldd,b];
Nout[l_,l1_,l2_,ldd_,b_]:=-(l+3+qv)Cm0[l+1]X[l+1,l1,l2,ldd,b]+(l-2-qv)Cm0[l]X[l-1,l1,l2,ldd,b]-(ldd-qv)Cm0[ldd+1]X[l,l1,l2,ldd+1,b]+(ldd+1+qv)Cm0[ldd]X[l,l1,l2,ldd-1,b];
Nmid[l_,l1_,l2_,ldd_,b_]:=(l1-1-qv)Cm0[l1+1]X[l,l1+1,l2,ldd,b]-(l1+2+qv)Cm0[l1]X[l,l1-1,l2,ldd,b]+(l2+2+qv)Cm0[l2+1]X[l,l1,l2+1,ldd,b]-(l2-1-qv)Cm0[l2]X[l,l1,l2-1,ldd,b];
FD[l_,l1_,l2_,ldd_]:=(X[l,l1,l2,ldd,Tanh[eta0+de]]-X[l,l1,l2,ldd,Tanh[eta0-de]])/(2 de);
comb[l_,l1_,l2_,ldd_]:=(1/2)Nout[l,l1,l2,ldd,bb]+(1/2)Nmid[l,l1,l2,ldd,bb]-bb X[l,l1,l2,ldd,bb];
report["Q1  cross-slot balance: (outer l,ldd shifts) = (middle l1,l2 shifts)",Max@Table[Abs[Nout[l,l1,l2,ldd,bb]-Nmid[l,l1,l2,ldd,bb]],{l,0,2},{l1,0,2},{l2,0,2},{ldd,0,2}],10^-22];
report["QA  N3:  Nout - beta X = d/deta X",Max@Table[Abs[FD[l,l1,l2,ldd]-(Nout[l,l1,l2,ldd,bb]-bb X[l,l1,l2,ldd,bb])],{l,0,1},{l1,0,1},{l2,0,1},{ldd,0,1}]];
report["QB  G5b: Nmid - beta X = d/deta X",Max@Table[Abs[FD[l,l1,l2,ldd]-(Nmid[l,l1,l2,ldd,bb]-bb X[l,l1,l2,ldd,bb])],{l,0,1},{l1,0,1},{l2,0,1},{ldd,0,1}]];
report["Q2  combined ODE (all three indices) vs d/deta X",Max@Table[Abs[FD[l,l1,l2,ldd]-comb[l,l1,l2,ldd]],{l,0,1},{l1,0,1},{l2,0,1},{ldd,0,1}]];
report["Q3  combined ODE on physical D (diagonal l1=l2=l')",Max@Table[Abs[FD[l,lp,lp,ldd]-comb[l,lp,lp,ldd]],{l,0,1},{lp,0,1},{ldd,0,1}]];
Print["--- all-index ODE harness done ---"];
