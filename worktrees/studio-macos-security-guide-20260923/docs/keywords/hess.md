# `[hess]`

The `[hess]` section controls Hessian, frequency, and thermochemistry
workflows.

## Keywords

### `type`

| Field | Value |
| --- | --- |
| Type | string |
| Default | `numerical` |
| Values | `numerical`, `analytical` |
| Used by | Hessian dispatch |

Selects numerical finite-difference or native analytical Hessian dispatch.
Analytical Hessians are validated for supported HF/DFT ground-state cases. The
input checker requests `type=numerical` when the method, state, or basis is
outside the analytical path.

### `state`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `0` |
| Used by | Hessian target state |

HF/DFT Hessians use state `0`. TDHF-family Hessian workflows use positive
excited-state indices.

### `dx`

| Field | Value |
| --- | --- |
| Type | float |
| Default | `0.01` |
| Used by | numerical Hessian displacement |

Finite-difference step size for numerical Hessians.

### `nproc`

| Field | Value |
| --- | --- |
| Type | integer |
| Default | `1` |
| Used by | numerical Hessian workers |

Worker count for numerical Hessian calculations. Must be at least `1`.

### `read`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | Hessian restart/input |

Reads an existing `.hess.json` sidecar where the workflow supports it. Current
sidecars carry a versioned state/model-configuration signature and the atom,
geometry, and isotopic-mass identity. OpenQP rejects unsigned older files or mismatches rather
than silently reusing a Hessian or vibrational data from another calculation;
rerun with `read=False` to generate a current sidecar.

### `restart`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | numerical Hessian restart |

Continues an interrupted numerical Hessian workflow.

Native TS/IRC drivers reject `restart=True` while their displaced-gradient
artifacts lack geometry/model signatures. Use it only for the standalone
numerical Hessian workflow; native TS/IRC initial Hessians are recomputed.

### `temperature`

| Field | Value |
| --- | --- |
| Type | float list |
| Default | `298.15` |
| Used by | thermochemistry |

Temperature or temperatures, in Kelvin, for thermochemical analysis. Values must
be positive.

### `clean`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | temporary Hessian files |

Removes temporary Hessian files where supported.

### `symmetry_unique`

| Field | Value |
| --- | --- |
| Type | boolean |
| Default | `False` |
| Used by | numerical Hessian displacements |

When enabled for a symmetric molecule, displaces one atom from each
symmetry-equivalent atom orbit and reconstructs the remaining numerical
Hessian rows with the detected abelian symmetry operations. OpenQP detects the
point group again at the current reference geometry and proves atom-orbit
coverage before launching displaced-gradient jobs.

If symmetry is absent, detection fails, the tolerance is too loose relative to
`dx`, or complete orbit coverage cannot be established, OpenQP reports the
reason and safely uses the full displacement set. The option does not enable
integral symmetry reduction in the displaced child calculations.
