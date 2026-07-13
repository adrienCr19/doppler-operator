(* ::Package:: *)
(* 25_tables_match.wl — do the CLOSED FORMS (bracket route, N8) reproduce the PREVIOUSLY PUBLISHED tables?
   Symbolic p-expansion of each closed form, substitute q -> -Ohat, compare to:
     (A) the site's Section-4 tables (the Ohat-polynomials displayed on index.html, from the ODE generator);
     (B) the literature SZ operator CR26 Eq. 11 (and its p^6 continuation).
   Fully symbolic; every check must reduce to 0. *)

g = Sqrt[1 + pp^2];                                   (* gamma in terms of p *)
Js[w_] := ((g + pp)^w - (g - pp)^w)/w;                (* symbolic J_w(q) *)
(* bracket table (verified numerically in harness 24), with a q-shift s: J_{n+q} -> J_{n+q+s} *)
F[0, 0, s_] := Js[2 + q + s]/(2 pp);
F[0, 1, s_] := Sqrt[3] (Js[1 + q + s] - g Js[2 + q + s])/(2 pp^2);
F[1, 0, s_] := Sqrt[3] (g Js[2 + q + s] - Js[3 + q + s])/(2 pp^2);
F[0, 2, s_] := Sqrt[5] (3 Js[q + s] - 6 g Js[1 + q + s] + (3 + 2 pp^2) Js[2 + q + s])/(4 pp^3);
F[2, 0, s_] := Sqrt[5] ((3 + 2 pp^2) Js[2 + q + s] - 6 g Js[3 + q + s] + 3 Js[4 + q + s])/(4 pp^3);
F[1, 1, s_] := 3 (g Js[1 + q + s] - (2 + pp^2) Js[2 + q + s] + g Js[3 + q + s])/(2 pp^3);
F[1, 2, s_] := Sqrt[15] (3 g Js[q + s] - 3 (3 + 2 pp^2) Js[1 + q + s] + g (9 + 2 pp^2) Js[2 + q + s] - (3 + 2 pp^2) Js[3 + q + s])/(4 pp^4);
F[2, 1, s_] := Sqrt[15] ((3 + 2 pp^2) Js[1 + q + s] - g (9 + 2 pp^2) Js[2 + q + s] + 3 (3 + 2 pp^2) Js[3 + q + s] - 3 g Js[4 + q + s])/(4 pp^4);
F[2, 2, s_] := 5 (3 (3 + 2 pp^2) Js[q + s] - 12 g (3 + pp^2) Js[1 + q + s] + 2 (27 + 24 pp^2 + 2 pp^4) Js[2 + q + s] - 12 g (3 + pp^2) Js[3 + q + s] + 3 (3 + 2 pp^2) Js[4 + q + s])/(8 pp^5);

(* closed form of D_{l lp ldd}, expanded in p to order ord, in Ohat = -q *)
Dcf[l_, lp_, ldd_, ord_] := Collect[Normal@Series[((-1)^(lp + ldd)/g) F[l, lp, 0] F[lp, ldd, -1], {pp, 0, ord}], pp, Simplify] /. q -> -OO;

chk[name_, expr_] := Print[name, ": ", If[Simplify[expr] === 0 || Simplify[expr] == 0, "PASS (== published)", "MISMATCH: " <> ToString[Simplify[expr]]]];

(* ---------- (A) the site's published Section-4 tables (lp = 0) ---------- *)
chk["A000", Dcf[0, 0, 0, 4] - (1 + (1/3) (OO - 3) OO pp^2 + (2/45) OO (OO^3 - 6 OO^2 + 5 OO + 12) pp^4)];
chk["A001", Dcf[0, 0, 1, 4] - (-((OO - 2) pp)/Sqrt[3] - (8 OO^3 - 42 OO^2 + 43 OO + 18) pp^3/(30 Sqrt[3]))];
chk["A002", Dcf[0, 0, 2, 4] - ((OO^2 - 5 OO + 6) pp^2/(3 Sqrt[5]) + (5 OO^4 - 42 OO^3 + 103 OO^2 - 42 OO - 72) pp^4/(63 Sqrt[5]))];
chk["A100", Dcf[1, 0, 0, 4] - ((OO - 1) pp/Sqrt[3] + (8 OO^3 - 30 OO^2 + 7 OO + 15) pp^3/(30 Sqrt[3]))];
chk["A101", Dcf[1, 0, 1, 4] - (-(1/3) (OO^2 - 3 OO + 2) pp^2 - (1/15) (OO^4 - 6 OO^3 + 7 OO^2 + 6 OO - 8) pp^4)];
chk["A102", Dcf[1, 0, 2, 4] - ((OO^3 - 6 OO^2 + 11 OO - 6) pp^3/(3 Sqrt[15]))];
chk["A200", Dcf[2, 0, 0, 4] - ((OO - 1) OO pp^2/(3 Sqrt[5]) + OO (5 OO^3 - 18 OO^2 - 5 OO + 18) pp^4/(63 Sqrt[5]))];
chk["A201", Dcf[2, 0, 1, 4] - (-OO (OO^2 - 3 OO + 2) pp^3/(3 Sqrt[15]))];
chk["A202", Dcf[2, 0, 2, 4] - ((1/45) OO (OO^3 - 6 OO^2 + 11 OO - 6) pp^4)];

(* ---------- (A) the site's published Section-4 tables (lp = 2, sample) ---------- *)
chk["A020", Dcf[0, 2, 0, 4] - ((1/45) OO (OO^3 - 6 OO^2 + 5 OO + 12) pp^4)];
chk["A222", Dcf[2, 2, 2, 4] - (1 + (1/21) (11 OO^2 - 33 OO - 54) pp^2 + (2/441) (23 OO^4 - 138 OO^3 - 80 OO^2 + 861 OO + 864) pp^4)];

(* ---------- (B) the literature SZ operator, CR26 Eq. 11 (and its p^6 term) ---------- *)
Sth = Collect[Normal@Series[((-1)^0/g) F[0, 0, 0] F[0, 0, -1] + (1/10) ((-1)^(2 + 0)/g) F[0, 2, 0] F[2, 0, -1] - 1, {pp, 0, 6}], pp, Simplify] /. q -> -OO;
chk["B-Sth p^2 (CR26 Eq.11)", Coefficient[Sth, pp, 2] - (OO^2 - 3 OO)/3];
chk["B-Sth p^4 (CR26 Eq.11)", Coefficient[Sth, pp, 4] - (7/150) OO (OO^3 - 6 OO^2 + 5 OO + 12)];
chk["B-Sth p^6 (extends CR26)", Coefficient[Sth, pp, 6] - (11/3150) OO (OO^5 - 9 OO^4 + 13 OO^3 + 57 OO^2 - 86 OO - 120)];

Print["--- tables-match harness done ---"];
