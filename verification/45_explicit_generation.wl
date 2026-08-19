(* ::Package:: *)
(* 45_explicit_generation.wl — FULLY EXPLICIT generation, no auxiliary coefficient families.
   Generates any D_{l l' l''} from the seed using ONLY:
       gamma = cosh(eta),  p = sinh(eta),
       J_w   = 2 sinh(w eta)/w                      (elementary bracket)
       C_l   = l/Sqrt(4 l^2 - 1)                    (explicit constants: 1/Sqrt3, 2/Sqrt15, ...)
       the seed D_000(q) = J_{q+1} J_{q+2} / (4 gamma p^2)
   and the four relations below.  NO U, rho, script-J, A, or any other named coefficient family
   is used or needed anywhere.  Every quantity on a right-hand side is either an explicit
   elementary scalar or a Doppler operator already generated at a lower index.

     (P0)  axis:      C_a D_{a00}(q) = (gamma/p) D_{a-1,00}(q)
                                       - (1/p)(J_{q+1}/J_{q+2}) D_{a-1,00}(q+1)
                                       - C_{a-1} D_{a-2,00}(q)
     (P3)  middle:    C_l' D_{h,l',0}(q) = (D_{l',00}(q-1)/p) [ gamma D_{h,l'-1,0}(q)/D_{l'-1,00}(q-1)
                                             - (J_{q-1}/J_q) D_{h,l'-1,0}(q-1)/D_{l'-1,00}(q-2) ]
                                           - C_{l'-1} (D_{l',00}(q-1)/D_{l'-2,00}(q-1)) D_{h,l'-2,0}(q)
     (N1)  reflection: D_{l l' l''}(q) = D_{l'' l' l}(-3-q)
     (P2)  fill:       D_{l l' l''}(q) = D_{l l' 0}(q) D_{l'' l' 0}(-3-q) / D_{0 l' 0}(q)

   The aberration kernel appears ONLY as the independent numerical check. *)

bb=42/100; gam=1/Sqrt[1-bb^2]; pp=gam bb; eta=ArcSinh[pp];
JJ[w_]:=2 Sinh[w eta]/w;                              (* elementary bracket *)
CC[l_]:=If[l<=0,0,l/Sqrt[4 l^2-1]];                   (* explicit constants *)
seed[q_]:=JJ[q+1] JJ[q+2]/(4 gam pp^2);               (* the only input *)

(* ---- (P0) axis tower, generated from the seed ---- *)
ax[0,q_]:=seed[q];  ax[a_,q_]/;a<0:=0;
ax[a_,q_]:=ax[a,q]=(1/CC[a])((gam/pp)ax[a-1,q]-(1/pp)(JJ[q+1]/JJ[q+2])ax[a-1,q+1]-CC[a-1]ax[a-2,q]);

(* ---- (P3) border tower, generated from the axis tower ---- *)
bd[h_,0,q_]:=ax[h,q];  bd[h_,lp_,q_]/;lp<0:=0;
bd[h_,lp_,q_]:=bd[h,lp,q]=(1/CC[lp])(
   (ax[lp,q-1]/pp)(gam bd[h,lp-1,q]/ax[lp-1,q-1]-(JJ[q-1]/JJ[q]) bd[h,lp-1,q-1]/ax[lp-1,q-2])
   -If[lp>=2,CC[lp-1](ax[lp,q-1]/ax[lp-2,q-1])bd[h,lp-2,q],0]);

(* ---- (P2)+(N1) interior fill ---- *)
gen[l_,lp_,ldd_,q_]:=bd[l,lp,q] bd[ldd,lp,-3-q]/bd[0,lp,q];

(* ---- independent check: the aberration-kernel definition ---- *)
Ybar[l_,x_]:=Sqrt[(2 l+1)/(4 Pi)]LegendreP[l,x];
Kk[d_?NumericQ,l_,lp_,s_]:=Kk[d,l,lp,s]=2 Pi NIntegrate[
   Ybar[l,mu]Ybar[lp,(mu+s bb)/(1+s bb mu)]/(gam(1+s bb mu))^d,{mu,-1,1},
   WorkingPrecision->30,PrecisionGoal->18];
ker[l_,lp_,ldd_,q_]:=(1/gam)Kk[-1-q,l,lp,-1]Kk[-q,lp,ldd,1];

qv=37/100;
Print["T1  explicit constants used: C_1..C_5 = ",Table[CC[l],{l,1,5}]];
Print["T2  axis tower from the seed (P0), a=0..5 vs kernel: ",
  ToString[N@Max@Table[Abs[ax[a,qv]-ker[a,0,0,qv]],{a,0,5}],InputForm]];
Print["T3  border tower from the axis tower (P3), h,l'<=5 vs kernel: ",
  ToString[N@Max@Table[Abs[bd[h,lp,qv]-ker[h,lp,0,qv]],{h,0,5},{lp,0,5}],InputForm]];
Print["T4  D_555 generated from the seed, vs kernel: ",
  ToString[N[Abs[gen[5,5,5,qv]-ker[5,5,5,qv]]],InputForm]];
Print["T5  ALL 216 elements generated from the seed, vs kernel: ",
  ToString[N@Max@Table[Abs[gen[l,lp,ldd,qv]-ker[l,lp,ldd,qv]],{l,0,5},{lp,0,5},{ldd,0,5}],InputForm]];
Print["T6  value check  D_555(q=0.37) = ",ToString[N[gen[5,5,5,qv],16],InputForm]];
Print["--- explicit-generation harness done ---"];
