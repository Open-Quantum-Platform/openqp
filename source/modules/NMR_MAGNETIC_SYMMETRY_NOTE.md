# Why the magnetic density response is exchange-only: J(P^B)=0 and f_xc[P^B]=0

**Scope.** Theoretical justification for the Phase-0 (ground-state) NMR coupled
magnetic response, and the foundation the MRSF-NMR extension rests on. It proves
that the first-order density response to a magnetic perturbation contributes
nothing through the Coulomb operator or through a semi-local exchange–correlation
kernel, so the only surviving two-electron coupling is exact (Hartree–Fock)
exchange. The accompanying numerical evidence is from the Phase-0 gates 0–2.

## 1. The first-order magnetic density is imaginary and antisymmetric

The external field `B` and the nuclear moment `m_N` enter through the magnetic
vector potential, giving purely **imaginary** one-electron perturbation operators
(orbital Zeeman `½L_O` and PSO `L_N/r_N³`, with `L = −i r×∇`). In a real AO basis
their matrices are real and **antisymmetric**. Standard first-order perturbation
theory then makes the first-order density-matrix response

```
    P^B  =  i · D ,      with D real and antisymmetric:  D^T = −D .
```

Equivalently `P^B` is anti-Hermitian; we carry the real antisymmetric matrix `D`
(the code's `P^B`) and verify `D^T = −D` numerically (gate 0).

## 2. The first-order charge density vanishes: ρ^B(r) = 0

The charge density built from any AO density matrix `M` is

```
    ρ_M(r) = Σ_{λσ} M_{λσ} φ_λ(r) φ_σ(r) .
```

The product `φ_λ(r) φ_σ(r)` is **symmetric** under `λ ↔ σ`. Contracting a symmetric
quantity with the antisymmetric `D` gives zero pointwise:

```
    ρ^B(r) = Σ_{λσ} D_{λσ} φ_λ(r) φ_σ(r) = 0      (for all r),
```

because `Σ_{λσ} (sym)_{λσ} (antisym)_{λσ} = 0`. **The first-order magnetic
perturbation produces no change in the (real) charge density.** This single fact
drives both vanishing-response results below.

## 3. Coulomb response vanishes: J(P^B) = 0

The Coulomb matrix is

```
    J_{μν}(M) = Σ_{λσ} (μν|λσ) M_{λσ} .
```

The two-electron integral `(μν|λσ) = ∫∫ φ_μφ_ν(1) r_12^{-1} φ_λφ_σ(2)` is
**symmetric under `λ ↔ σ`** (the second charge distribution `φ_λφ_σ` is symmetric
and `r_12^{-1}` is symmetric). Contracting with antisymmetric `D`:

```
    J_{μν}(D) = Σ_{λσ} (μν|λσ) D_{λσ} = 0      (every μν).
```

Physically: `J` is the classical Coulomb potential of `ρ^B`, and `ρ^B = 0`
(Section 2), so there is no first-order Coulomb potential. **Proven, element-wise.**

## 4. Semi-local XC-kernel response vanishes: f_xc[P^B] = 0

A semi-local (LDA/GGA/meta-GGA) functional makes the XC potential a function of
the local, real density variables — `ρ`, `∇ρ`, and the kinetic-energy density `τ`.
Its linear response to a density-matrix perturbation `M` is

```
    f_xc[M](r) = ∫ K_xc(r, r') δρ_M(r') dr'      (+ gradient/τ terms),
```

i.e. it is driven entirely by the **density change** `δρ_M`, `∇δρ_M`, `δτ_M`. For
`M = P^B` we showed `ρ^B(r) ≡ 0`; the same symmetry argument applied to
`∇(φ_λφ_σ)` and to `Σ ∇φ_λ·∇φ_σ` (both symmetric in `λ↔σ`) gives `∇ρ^B = 0` and
`δτ^B = 0`. Hence every semi-local response channel vanishes:

```
    f_xc[P^B] = 0      for any ρ/∇ρ/τ-dependent functional.
```

**Caveat (why this is exactly the standard DFT-NMR approximation).** The magnetic
perturbation *does* induce a first-order paramagnetic **current** `j_p`. A
current-density functional `E_xc[ρ, j_p]` would respond to it. Standard functionals
omit the current dependence, so `f_xc[P^B] = 0` is exact *within* that
(conventional) class — and any future current-DFT term would enter precisely here.

## 5. Exact-exchange response survives: K(P^B) ≠ 0

The exchange matrix is

```
    K_{μν}(M) = Σ_{λσ} (μλ|νσ) M_{λσ} .
```

Here `λ` is paired with `μ` and `σ` with `ν`, so `(μλ|νσ)` has **no `λ↔σ`
symmetry**; the contraction with antisymmetric `D` does not vanish. Moreover
`K(D)` is itself antisymmetric, preserving the structure of the response. This is
the lone two-electron coupling in the magnetic CPHF/CPKS equations, and it is
scaled by the functional's exact-exchange fraction `c_x` (1 for HF, 0 for pure
GGA, 0.5/0.25 for BHHLYP/PBE0). Hence for pure functionals the coupled response
collapses to the uncoupled one, while HF/hybrids retain an exchange coupling
proportional to `c_x`.

## 6. Numerical evidence (Phase-0 gates 0–2)

H2O / STO-3G, common gauge origin at the center of mass; `P^B` is the converged
first-order density (z field component); norms are Frobenius.

| Quantity | HF | BHHLYP | PBE0 | PBE | expected |
|----------|----|--------|------|-----|----------|
| gate 0 `max\|P^B + (P^B)^T\|` | 1.1e-16 | 1.1e-16 | 1.1e-16 | 1.1e-16 | 0 (antisymmetric) |
| gate 1 `\|\|J(P^B)\|\|`        | 6.6e-17 | 1.2e-16 | 2.9e-16 | 1.5e-16 | 0 (Coulomb vanishes) |
| gate 2 `\|\|K(P^B)\|\|` (c_x=1)| 0.5417  | 0.5417 | 0.5417 | 0.5417 | > 0 (exchange survives) |

End-to-end corollary (gate 3): for **pure PBE** the coupled and uncoupled
shieldings are *identical to machine precision*, which is only possible if both
`J(P^B)` and `f_xc[P^B]` contribute nothing — a direct experimental confirmation
of Sections 3–4. Conversely (gates 4, 6) the HF/hybrid coupled − uncoupled
difference is nonzero and scales monotonically with `c_x`, confirming Section 5.

## 7. Consequence for MRSF-NMR

The same symmetry holds for the MRSF magnetic response: the first-order density is
still imaginary/antisymmetric, so the MRSF coupled response inherits "Coulomb-zero,
semi-local-fxc-zero, exact-exchange-only." The MRSF-specific content is therefore
confined to (i) the spin-flip/mixed-reference structure of the relaxed density and
(ii) the amplitude/Z-vector relaxation — not to new two-electron response channels.
This is the foundation for the gate-5a ground-state-limit check in
`MRSF_NMR_DESIGN.md`.
