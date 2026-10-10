# QMRSF-DK

QMRSF-DK builds a spin-adapted response space from a **quintet** reference and
dresses its closed-shell block with exact exchange. From one quintet ROHF/ROKS
calculation it returns the complete active space of four electrons in four
orbitals, CAS(4,4): 20 singlets, 15 triplets, and one quintet.

## Background

The four singly occupied orbitals of a quintet reference define the active
space. Single spin flips from the eight complementary `Ms = ±1` reference
determinants generate every two- and four-open-shell configuration; the two
`Ms = ±2` determinants supply the six closed-shell configurations that single
spin flips cannot reach, through a dressed kernel. The resulting blocks are
small enough to diagonalize directly, so no Davidson iteration is needed and
the response costs little beyond the reference SCF.

The method targets systems whose frontier manifold is strongly correlated:
open-shell singlet ground states, antiaromatic transition states, and the
triplet-pair states of singlet fission, where the singlet, triplet, and quintet
descriptions of the same pair must be available on an equal footing.

## Exchange partition and the S3R default

A hybrid functional assigns two-electron integrals to the exact-exchange class
by their orbital indices. That assignment is not preserved by rotations inside
a degenerate active-orbital pair, so a reference with an exactly or nearly
degenerate frontier pair gives excitation energies that depend on the
orientation the SCF happened to converge to.

QMRSF-DK evaluates three exchange partitions in the same run:

| Convention | Description |
| --- | --- |
| `value-based` | the original index-based assignment; orientation dependent |
| `S3R` | the covariant convention used for the reported result |
| `Haar` | the average over the full pair-rotation group; also covariant |

**S3R is the production result.** The value-based and Haar spectra are printed
alongside it as diagnostics, together with a covariance gate: the largest
change of any excitation energy when the detected degenerate pairs are rotated.
The covariant conventions return machine zero for that gate; the value-based
one exposes the orientation freedom, which reaches several electronvolts in
degenerate systems.

Pairs are detected from the SOMO energies. Every adjacent pair degenerate to
within `1e-2` hartree is averaged, which covers symmetry-enforced degeneracies
and the distance-degenerate frontier pairs of a separated dimer. A reference
with no such degeneracy has no orientation freedom to remove, and the same
construction is then applied to the frontier pair alone, so the definition
stays continuous as a degeneracy is approached.

## Requirements

| Setting | Required value |
| --- | --- |
| `[input] runtype` | `energy` |
| `[input] method` | `tdhf` |
| `[scf] type` | `rohf` |
| `[scf] multiplicity` | `5` |
| `[tdhf] type` | `qmrsf_dk` |
| `[tdhf] multiplicity` | `1` |

The pre-flight checker reports any of these that is missing before the
calculation starts.

## Minimal example

Square H4, whose degenerate frontier pair is exactly what the covariant
partition is for:

```ini
[input]
system=
   1   0.700000000   0.700000000   0.000000000
   1  -0.700000000   0.700000000   0.000000000
   1  -0.700000000  -0.700000000   0.000000000
   1   0.700000000  -0.700000000   0.000000000
charge=0
runtype=energy
basis=sto-3g
functional=bhhlyp
method=tdhf

[guess]
type=huckel

[scf]
type=rohf
multiplicity=5

[tdhf]
type=qmrsf_dk
hfscale=0.5
```

`examples/QMRSF-DK/H4_QMRSF-DK_S3R_ENERGY.inp` in the code repository is this
input.

## Range-separated references

A range-separated functional treats `alpha + beta*erf(mu*r12)` of the electron
interaction by exact exchange, so its exchange operator is

```
K = alpha * K[1/r12] + beta * K[erf(mu*r12)/r12]
```

rather than a single scaled exchange. QMRSF-DK uses that operator in both
places where the exchange fraction enters, the reference core exchange and the
dressed kernel, and applies the covariant partition to the whole operator
rather than to its short-range part alone. Setting `beta = 0` reproduces the
global-hybrid result exactly.

```ini
[input]
functional=camb3lyp
...

[tdhf]
type=qmrsf_dk
```

The reference `alpha`, `beta`, and `mu` come from the functional. `[tdhf]
cam_alpha`, `cam_beta`, and `cam_mu` override them for the kernel; a negative
value inherits the corresponding reference value independently. Only a
resolved set that remains incomplete or invalid (for example, a nonpositive
range parameter when range separation is active) is refused.

Orientation dependence is *larger* with a range-separated reference than with a
global hybrid, because more of the exchange operator carries the noncovariant
assignment, so the covariant partition matters more there.

## Output

The log reports the three spectra, the leading configurations of each state,
the covariance gate, and the detected degenerate pairs with their doublet
anisotropy. Machine-readable results are written next to the log as
`<log>.qmrsf_dk.json`, which carries the S3R singlet, triplet, and quintet
states with their total energies and excitation energies. Every manifold is
reported against the lowest singlet, so the singlet-triplet and
singlet-quintet gaps appear directly. The singlet excitation energies are also
published through the standard `OQP::td_energies` result.

`qmrsf_dk_active.dat` holds the active-space integrals for an independent
check of the same construction. It carries the plain integrals only, so it does
not reproduce a range-separated spectrum on its own.

## Limitations

- QMRSF-DK is currently invoked through sectioned `.inp` input. The concise
  `.oqp` grammar and the high-level `OpenQP` Python helpers do not yet expose
  `[tdhf] type=qmrsf_dk`.
- Energies only; gradients and other runtypes are not implemented.
- The response covers the active `O -> O` space. Excitations that need
  correlation from outside the four frontier orbitals, such as the ionic
  `pi pi*` states of benzene, are overestimated.
- The active space is fixed at CAS(4,4) by the quintet reference.
