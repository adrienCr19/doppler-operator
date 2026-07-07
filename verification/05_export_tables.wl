Cm[l_, m_] := If[l <= Abs[m], 0, Sqrt[(l^2 - m^2)/(4 l^2 - 1)]];
buildC[lp_, m_] := Module[{c},
  c[k_][l_, ldd_] /; l < Abs[m] || ldd < Abs[m] := 0;
  c[0][l_, ldd_] := If[l == lp && ldd == lp, 1, 0];
  c[k_][l_, ldd_] := c[k][l, ldd] = Expand[(1/k) (
      -(l + 3 + q) Cm[l + 1, m] c[k - 1][l + 1, ldd] + (l - 2 - q) Cm[l, m] c[k - 1][l - 1, ldd]
      - (ldd - q) Cm[ldd + 1, m] c[k - 1][l, ldd + 1] + (ldd + 1 + q) Cm[ldd, m] c[k - 1][l, ldd - 1]
      - Sum[SeriesCoefficient[Tanh[x], {x, 0, j}] c[k - 1 - j][l, ldd], {j, 1, k - 1}])];
  c];
toP[c_, l_, ldd_, ordE_] := Collect[Normal@Series[Sum[c[k][l, ldd] et^k, {k, 0, ordE}] /. et -> ArcSinh[pp], {pp, 0, ordE}] /. q -> -OO, pp, Simplify];
c04 = buildC[0, 0]; c24 = buildC[2, 0];
fmt[e_] := StringReplace[ToString[TeXForm[e /. OO -> Global`\[ScriptCapitalO]]], "\\text{$\\scriptcapitalo$}" -> "\\hat{O}"];
fmt2[e_] := StringReplace[ToString[TeXForm[e]], {"\\text{OO}" -> "\\hat{O}_\\nu", "OO" -> "\\hat{O}_\\nu", "\\text{pp}" -> "p", "pp" -> "p"}];
tab = Join[
   Flatten[#, 1] &@Table[{l, 0, ldd, fmt2[toP[c04, l, ldd, 4]]}, {l, 0, 2}, {ldd, 0, 2}],
   Flatten[#, 1] &@Table[{l, 2, ldd, fmt2[toP[c24, l, ldd, 4]]}, {l, 0, 4}, {ldd, 0, 4}]];
Export["doppler_tables_O.json", Map[<|"l" -> #[[1]], "lp" -> #[[2]], "ldd" -> #[[3]], "tex" -> #[[4]]|> &, tab], "JSON"];
(* S_th assembly to p^6 *)
c06 = buildC[0, 0]; c26 = buildC[2, 0];
toP6[c_, l_, ldd_] := Collect[Normal@Series[Sum[c[k][l, ldd] et^k, {k, 0, 6}] /. et -> ArcSinh[pp], {pp, 0, 6}] /. q -> -OO, pp, Simplify];
S = Collect[Expand[toP6[c06, 0, 0] + toP6[c26, 0, 0]/10 - 1], pp, Simplify];
Print["Sth p2: ", fmt2[Coefficient[S, pp, 2]]];
Print["Sth p4: ", fmt2[Coefficient[S, pp, 4]]];
Print["Sth p6: ", fmt2[Coefficient[S, pp, 6]]];
Print["done"];
