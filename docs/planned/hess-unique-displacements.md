# Symmetry-unique displacements in the numerical Hessian

## What

`pyoqp/oqp/library/single_point.py:1385` builds
`shift = np.diag(np.ones(ncoord)*dx)` — all 3N coordinates, both directions,
i.e. 6N gradient jobs. Symmetry-related displacements are exactly redundant:

    H[σ(a)α, σ(b)β] = Σ M_αγ M_βδ H[aγ, bδ]

This is **exact, not an approximation** — a symmetry-breaking imaginary mode is
still recovered in full.

## Orbit reduction

| molecule | atoms → orbits | factor |
|---|---|---|
| H2O | 3 → 2 | 1.5x |
| C2H4 | | 3x |
| naphthalene | 18 → 5 | 3.6x |
| benzene | 12 → 2 | 6x |

**But only with the full non-abelian group** from `enumerate_full_group`; the
abelian subgroup staged today gives benzene 4 orbits, i.e. 3x.

## Two honest deflations

1. **Wall clock is capped by the pool width.** `ncpu = np.amin([ndim, nproc])` at
   `:1414`, so the real gain is `ceil(6N/ncpu_full) / ceil(6N_red/ncpu_red)` —
   water on ≥18 cores gains nothing. Realistic target: 2–4x on mid-size
   symmetric molecules.
2. **The extra ± halving needs the standard orientation**, which `numerical_hess`
   is explicitly excluded from (`molecule.py:445-447`). In a random input frame
   it delivers exactly zero. Recovering it needs site-group-adapted displacement
   directions plus a 3×3 back-transform.

## The unconditional win

The **serial** IR/Raman dipole-derivative loop (`:1202-1216`), where the CPU
factor *is* the wall-clock factor and each displaced evaluation is a full
property run containing a 3-RHS CPHF.

## Verification

Reduced vs full Hessian to < 1e-8 Ha/Bohr² and frequencies to < 0.1 cm⁻¹ on
H2O/C2v, C2H4/D2h, benzene/D6h; one deliberately symmetry-broken case must fall
back to 6N. Assert the geometry possesses the group to much better than `hess.dx`
(default 0.01 Bohr) — refuse if `tolerance > dx/10`. Two symmetry-related
displacements must converge to symmetry-related SCF solutions (not guaranteed for
ROKS or when stability following kicks in); for MRSF the state irrep label must
agree across an orbit. Keep `hess.restart` working — the on-disk gradient indices
change meaning.

## Ordering

**The child-reorientation fix must land first**, or this builds on a wrong
Hessian.
