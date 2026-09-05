# State-overlap minor determinants for NACME and NAMD (`tdhf(tlf=...)`)

The MRSF state overlap between consecutive geometries,
`<Psi_I(t-dt)|Psi_J(t)>`, is evaluated by determinant factorization: the
reference determinant is shared, so the overlap reduces to a contraction of the
response amplitudes with three classes of minor determinants of the MO overlap
matrix (`s_ij` one-hole occupied minors, `s_ab` particle minors, `s_ia` mixed
minors) rather than one Loewdin determinant per amplitude pair.

`tlf` selects how the `s_ij` and `s_ab` minors are evaluated (`s_ia` is always
exact):

| `tlf` | minors | notes |
| --- | --- | --- |
| `0` = `notlf` = `exact` (default) | exact Gaussian-elimination minors, no truncation | independent of orbital rotations between steps; this is *not* the paper's zeroth-order TLF(0), which is not implemented |
| `1` | first-order truncated Leibniz formula, TLF(1) | JCTC 15, 882 (2019) |
| `2` | second-order truncated Leibniz formula, TLF(2) | most accurate TLF approximation; KNU-GAMESS `ndtlf=2` |

The truncated Leibniz formula assumes that the MOs of consecutive steps are
nearly orthonormal, i.e. that the MO overlap matrix is close to diagonal.  When
near-degenerate doubly occupied orbitals rotate into each other within one
nuclear step (a 45-degree mixing of two occupied orbitals was observed in hot
uracil trajectories), the diagonal MO overlaps drop to ~0.7 and TLF(2) returns a
collapsed state overlap (all diagonal elements ~0.3-0.4) even though the SCF
solution and the MRSF surfaces are continuous.  Norm-preserving interpolation
then turns the collapsed overlap into a large spurious time-derivative coupling.
The exact minors are invariant to such rotations and, for molecules of the size
of uracil (30 occupied alpha orbitals, 6-31G*), cost the same wall time as
TLF(2).  For large systems where the `nvir^2` particle minors dominate, the
recommended route is Jacobi's complementary-minor identity (all one- and two-hole
minors from one LU factorization of the occupied block) rather than truncation.

The state-overlap section of the log states which evaluation was used:
`state-overlap minors: exact minor determinants (tlf=0, default; ...)` or
`TLF(n) truncated-Leibniz minors; ...`.  NAMD additionally logs a warning when
every column norm of the retained state overlap falls below 0.5.
