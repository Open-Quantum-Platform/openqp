# PCM/ddX

The current PCM production path is an energy-only reference-SCF workflow using
the ddX backend.

```ini
[input]
runtype=energy
method=hf
basis=6-31g*
ispher=true

[scf]
type=rhf
multiplicity=1

[pcm]
enabled=true
backend=ddx
mode=reference_scf
model=ddpcm
epsilon=78.3553
```

Runnable reference: `examples/PCM/H2O_RHF-HF_DDPCM_ENERGY_ISPHER.inp`.

## Scope

- Supported production mode: `backend=ddx`, `mode=reference_scf`,
  `runtype=energy`.
- Reference support: RHF and ROHF.
- MRSF-TDDFT can use the high-spin ROHF reference density as a reference-SCF
  PCM baseline.
- PCM gradients, PCM optimizations, Hessians, NACs, and state-specific
  excited-state PCM are outside this first energy path.
