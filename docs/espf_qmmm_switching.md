# ESPF QM/MM grid switching: `ESPF_SWSCALE` at covalent boundaries

The ESPF electrostatic-embedding gradient carries a small finite-difference
residual (~1e-3 a.u., the grid-derivative floor). At a **covalent** QM/MM
boundary that residual is not spread over the QM region: it concentrates as an
over-shoot on the MM host atom (M1), and it is controlled by the smooth-
switching parameter `ESPF_SWSCALE` — the grid weight-derivative term in
`espf_grad_weight` (`source/modules/qmmm.F90`).

The shipped default, `ESPF_SWSCALE=1.8`, **over-smooths at covalent
boundaries**.

## Recommendation

| QM/MM boundary | recommended `ESPF_SWSCALE` |
| --- | --- |
| covalent (a bond is cut; link atoms present) | **1.5** |
| whole-molecule (no bond cut) | 1.8 (the default) |

`ESPF_SWSCALE` is read from the environment, so adopting the recommendation
needs no rebuild:

```bash
export ESPF_SWSCALE=1.5
```

The default is deliberately left at 1.8: 1.5 is better at covalent boundaries
and worse without one, so it is not a safe global value.

## Measurement behind it

Finite-difference check of the analytic force against the numerical `-dE/dx` of
the *full self-consistent* QM/MM energy (`QMMM_MD.oqp_driver.compute_force`,
HF/6-31G, AMBER-14 MM, `frontier_scheme=none`). Worst `|analytic - FD|` over QM
atoms and link hosts, in kJ/mol/nm:

| System | Boundary | `SWSCALE=1.8` | `SWSCALE=1.5` |
| --- | --- | --- | --- |
| alanine dipeptide, QM = C-terminal amide (cut ALA C–CA) | link atom | 14.7 | **5.1** |
| alanine dipeptide, QM = ACE cap (cut ACE C–ALA N) | link atom | 67.8 | **13.0** |
| water dimer, QM = one water | none | **6.8** | 8.2 |

So 1.5 cuts the boundary-host residual by 5–13x on both covalent-boundary
systems, and is about 20% worse for the whole-molecule case.

A parallel sweep of `ESPF_SWDELTA` against `ESPF_SWSCALE` on the amide system
showed `ESPF_SWSCALE` is the dominant knob: 1.5 is far better than 1.8, while
2.1 diverges to ~1000s. `ESPF_WSCALE` must stay at 1.0.

## Related knobs

All are runtime environment variables read in `source/modules/qmmm.F90` and
mirrored in `pyoqp/oqp/library/namd.py`:

| Variable | Default | Meaning |
| --- | --- | --- |
| `ESPF_SMOOTH` | on | enable smooth grid switching at all |
| `ESPF_SWSCALE` | 1.8 | switching width scale (see above) |
| `ESPF_SWDELTA` | 0.7 | switching offset |
| `ESPF_WSCALE` | 1.0 | grid weight scale; keep at 1.0 |
| `ESPF_WDERIV` | on | include the grid weight-derivative gradient term |

## Scope

This is a property of baseline ESPF and therefore affects all QM/MM, including
SOC-NAMD-QMMM. It is independent of the frontier-charge work: with the RCD
frontier scheme the analytic gradient already matches the full-field baseline,
and is better at the boundary host because RCD removes the raw M1 charge.

Broader multi-system validation is welcome before the shipped default is
reconsidered; see issue #260.
