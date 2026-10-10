# OQP Studio Troubleshooting

## The runner is unavailable

In Execution, inspect the runner list, detected version, and executable path.
**OpenQP (local)** means an `openqp` executable found on Studio's search path;
it does not mean the engine bundled in the application. Select **OpenQP
(bundled)** for a `with-engine` installation.

If a pip-installed OpenQP is not found, compare the path shown by Studio with:

```bash
command -v openqp
openqp --version
```

Applications started from Finder or the Windows desktop do not necessarily
inherit the same shell `PATH`. Use the bundled runner or install OpenQP in a
standard executable location visible to desktop applications.

## A calculation does not start

Check, in order:

1. a runner is selected and shown as available;
2. the input filename ends in `.oqp` or `.inp`;
3. the results folder shown in Execution is writable;
4. the memory admission message permits the calculation; and
5. the Live log for the first OpenQP error.

On macOS, a denied Documents-folder prompt can make Studio choose a fallback
result directory. Select an explicit folder from **File > Results folder**.

## The calculation failed or did not converge

**failed** means the engine returned an error or could not be launched. Read the
project's calculation log in Project files; the Live log is also retained in
the project directory.

**not converged** means Studio recognized an incomplete optimization and kept
the best available geometry. Select **Restart** to prepare a new calculation
from that geometry. Increasing iterations or changing trust and convergence
controls should be a deliberate scientific choice based on the log, not an
automatic response to every failure.

For a transition state, completion of an optimizer is not sufficient evidence
of a transition state. Calculate the Hessian, verify one imaginary frequency,
and use IRC when the reaction connection must be established.

## PCM reports missing ddX

The message

```text
PCM (ddX) cavity build failed: OpenQP was built without OQP_ENABLE_DDX
```

comes from the selected OpenQP engine. Choose a bundled engine that includes
ddX on macOS or Linux, or rebuild a local engine with ddX enabled on a supported
platform. Native Windows OpenQP builds currently reject ddX; select a
ddX-enabled WSL runner for PCM on Windows. This is independent of the Studio
input form.

## An orbital or normal mode does not appear

SCF orbitals require a compatible Molden file and appear only after an orbital
is selected. Dyson orbitals require an IP/EA calculation and its dedicated
Dyson orbital output. NTOs additionally require saved state-resolved transition
density data. Use **Project files** to confirm that the expected JSON and
Molden files exist.

Normal modes require both frequencies and displacement vectors. Select a mode
from **Choose a normal mode**; select **Reset** if an earlier surface or
animation state obscures the equilibrium molecule.

## A spectrum is absent

The Spectrum card is data driven. Absorption needs excitation energies and
oscillator strengths. ESA needs transitions from the selected excited state to
higher states. Emission should be interpreted at an optimized excited-state
structure with a corresponding downward transition. IR needs frequencies and
intensities, while photoelectron spectra need EKT Dyson roots.

Open the calculation log and confirm that the relevant table was actually
written. A method name alone does not guarantee that every spectrum can be
constructed.

## macOS startup diagnostics

The desktop app normally starts its backend internally. To isolate a packaged
backend startup problem, run the sidecar from Terminal so its standard error is
visible:

```bash
time "/Applications/OQP Studio.app/Contents/MacOS/oqp-studio-backend" --port 8899
```

This diagnostic command uses the backend's optional HTTP mode; it does not
describe normal desktop operation. Startup timing lines beginning with
`startup` separate frozen-runtime extraction, Python imports, and service setup.
The failure-only backend log is normally under:

```text
~/Library/Application Support/OQP Studio/backend.log
```

When reporting a startup problem, include the complete standard error, timing
lines, macOS version, machine architecture, Studio asset name, and whether a
Documents access dialog appeared.

## Network lookups fail

PubChem, PDB, update checks, and on-demand engine installation require network
access. Coordinate editing and bundled-engine calculations do not. On a campus
or company network that re-signs TLS traffic, open **File > Network settings**
and add the institution's trusted PEM or CRT root certificate. Do not import a
certificate obtained from an unverified source.

For OpenQP engine errors unrelated to the desktop interface, continue with the
general [OpenQP troubleshooting guide](../troubleshooting.md).
