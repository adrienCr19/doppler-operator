# The Doppler Operator Project

A machine-verified research catalog for the **boost operator** of relativistic Compton-scattering theory
(Chluba & Ravenni 2026; Chluba & Rosenberg 2026; Rosenberg & Chluba 2026; Hoey, Long & Chluba 2026) and a
complete **generation calculus for the Doppler operator**
`D^m_{ℓℓ'ℓ''}(ν,β) = γ⁻¹ · ⁻¹B(ν,−β) · ⁰B(ν,β)` — the round-trip scattering combination whose coefficient
structure the literature explicitly left as "a problem … for future work" (Rosenberg & Chluba 2026, §2.3).

**Live site:** https://adriencr19.github.io/doppler-operator/ ·
**Manuscript:** [`manuscript/paper.pdf`](manuscript/paper.pdf) ·
**Release bundle:** [`doppler-operator-bundle.zip`](doppler-operator-bundle.zip)

## Results (N-catalog)

**Centerpiece — the exact recurrence relations** (`proofs/recurrences.html`): the index-raising relations of the
Doppler operator, each proved by descending to the aberration kernel and carrying every line of algebra to its
conclusion. The four kernel identities (seed `ᵈK₀₀ = J₁₋d/2p`, parity, transposition, and the master ℓ-raising
recurrence from the Legendre relation) are derived from the integral definition; N5/N6/(A) are reassembled from
them. All steps verified to deviation 0 (harness 36).

| | |
|---|---|
| N1/N1a, N2 | reflection `D(Ô) = Dᵀ(3−Ô)` (diagonals are polynomials in `D̂_ν = Ô²−3Ô`), parity |
| N3, N3′ | closed rapidity ODE (the table generator) and the all-index combined flow with cross-slot balance |
| N4–N7 | sum rule; outer-index raisings on the weight lattice; seeds + finite generation algorithm + `₂F₁` closed form |
| N8 | closed forms as bracket products of `J_w = 2 sinh(wη)/w` (matches all published tables symbolically) |
| N9–N11 | thermal function class (iterated Bessel antiderivatives), middle-index raising, the F-calculus |
| N12–N14 | reflection = self-adjointness on photon phase space; Henyey–Greenstein generating identity; recoil rungs = weight lattice |
| N15 | polarized sector foundations; RC26 Eq. 20 verified operator-exactly, shown collective (fails per-m) |
| N16 | kSZ dipole: `D₁₂₁` closed form, exact `M₁₁` tower, reflection **obstruction theorem** |
| N18 | `SL(2,ℂ)` dictionary: boosts = celestial dilations, unitarity on `ν^{1−2d}dν`, principal series at complex weight |
| N19 | **proof of RC26 Eq. 20**: exchange symmetry `ₛK^{d,m} = ₘK^{d,s}` + on-curve swap lemma (machine certificates) |
| N20 | Compton dispersion: marginal mode = boosted comoving equilibrium; **exact flatness** `T(k) = W_k T(0) W_k⁻¹`; thermal gap |
| N21 | frame-exact kinematic `C_k`: reproduces Nozawa's `C₁` symbolically, proves the `τ/τ*` shift `= +D̂_ν`, new `C₂, C₃` |
| N22/N23 | **rationalization lemma**: every kernel element (any spin, any m) is a finite J-combination; polarized catalog |
| N24 | polarized ladders: spin-weighted N5/N6 and the middle-index relation `M_s` on the mixed-spin family |
| N25 | **raising every index with Doppler operators only**: pure-`D` middle-index step (P3) via the factor dictionary, plus axis ladders (P0), reflection (N1), rank-one fill (P2) — no kernels, brackets, or extended `X` family |
| E1 | the exact thermal SZ operator assembled in closed form (elementary + Bessel + one transcendent with full calculus) |

## Contents

- `index.html` — the catalog (§1 boost operator with citations, §2 Doppler operator, §3 core recursions,
  §4 structure theory, §5 sectors & applications, §6 tables, §7 verification).
- `proofs/` — one self-contained proof page per result (see site navigation).
- `verification/` — 36 Wolfram-language harnesses + saved PASS logs (`out_*.txt`) + browsable per-script pages.
- `data/` — machine-readable exact tables:
  `doppler_tables_p8.json` (501 scalar elements, O(p⁸), exact q-polynomial coefficients) and
  `polarized_tables_p6.json` (81 polarized elements T(2,2), T(0,2), T(2,0), m = 0…2, O(p⁶)).
- `manuscript/` — the companion paper (LaTeX + PDF).

## Reproduce

```
cd verification
wolframscript -file 01_conventions.wl     # conventions + literature identities
wolframscript -file 02_doppler.wl         # N1–N6
...                                       # see verification/index.html for the full list
wolframscript -file 30_spin_m_closed_forms.wl
```

Exact identities pass at the 1e−20…1e−27 level (25–40-digit quadrature); ODE checks are finite-difference
limited (~1e−12); truncated-sum checks are tail-limited. Key results (closed forms vs. published tables,
kinematic C₁, obstruction/gap factorizations, on-curve swap certificates) are **symbolic zeros**.

## Status

Research notes: machine-verified, not peer-reviewed. Derivations and verification prepared with
Claude (Anthropic) + Wolfram 15.
