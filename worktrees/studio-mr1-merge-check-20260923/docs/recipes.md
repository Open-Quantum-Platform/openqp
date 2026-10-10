# End-to-end recipes

An OQP Studio recipe is a portable JSON document with schema
`oqp-studio-recipe/v1`. It applies a scientific workflow and preserves the
settings needed after input generation:

- `method`: OpenQP method, basis, state, and workflow-specific controls
- `execution`: runner preference and thread count
- `analysis`: spectrum, state, width, and property-map choices
- `art`: content, quality, exposure, background, depth of field, and surface appearance
- `postprocess`: optional Python run after a successful calculation

The molecule and result directory are deliberately excluded. A recipe can
therefore be applied to a new structure and stores no machine-specific path.
Unknown setting names are ignored so a newer recipe can still be opened by an
older Studio without altering unrelated controls.

```json
{
  "schema": "oqp-studio-recipe/v1",
  "id": "vertical-absorption",
  "name": "Vertical absorption and NTO",
  "description": "Calculate vertical excited states and prepare NTO analysis.",
  "tags": ["photochemistry", "absorption", "NTO"],
  "available": true,
  "unavailable_reason": "",
  "workflow": "abs",
  "method": {"theory": "mrsf", "nstate": 6},
  "execution": {"threads": 4},
  "analysis": {
    "spectrum_kind": "absorption",
    "spectrum_shape": "lorentzian",
    "spectrum_width": 20
  },
  "art": {
    "content": "excited:nto_hole",
    "quality": "0.75,5",
    "background": "studio"
  },
  "postprocess": {"enabled": false, "code": "", "timeout_seconds": 60}
}
```

## Python after calculation

Python embedded in an imported recipe is disabled until its source is opened
under **Review Python** and trusted on the current computer. Local trust is
keyed to the recipe ID and exact source hash. It is stored outside the recipe,
so exporting the same JSON cannot grant trust on another computer.

Approved code runs in a separate process with the calculation directory as
its current directory. `JOB_DIR` is a `pathlib.Path` for that directory.
Standard output and standard error are written to `postprocess.log`, and the
process is terminated at the recipe timeout. This is not a Python sandbox:
trusted code can read or modify files and use the network.

## NICS

NICS requires magnetic shielding at a non-nuclear probe center and is defined
as the negative isotropic shielding at that point. OpenQP currently emits NMR
shielding only for real nuclei, so the built-in NICS recipe remains visible but
cannot be applied. It can be enabled when OpenQP accepts ghost/probe centers
and includes their shielding values in the normal JSON result.
