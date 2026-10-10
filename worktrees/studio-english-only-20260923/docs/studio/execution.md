# Execution and Projects

## Review the input

Execution contains the generated OpenQP input and its filename. The filename
must end in `.oqp` or `.inp` and cannot contain a directory path. Studio normally
uses the project name for the input and log filenames.

The editor accepts an existing concise `.oqp` input and, for compatibility,
legacy `.inp` text. Workflow generation itself produces concise `.oqp` input.

## Select where calculations run

Choose a runner before **Run** becomes available:

- **OpenQP (bundled)** uses the engine delivered with or installed by Studio.
- **OpenQP (local)** uses an `openqp` executable found on the local search path.
- **WSL** uses a detected OpenQP installation in Windows Subsystem for Linux.

The menu displays each engine's detected version when available. Selecting a
runner is explicit; Studio does not silently substitute the bundled engine for
a missing local installation.

## Hardware and memory

**Execution host** reports logical and physical CPU information, total and
available memory, and the active platform. Set **OpenMP threads** according to
the calculation and available physical cores.

Studio estimates memory demand from the input before submission. If the
estimate exceeds currently available RAM, the job is refused with a RAM-limit
message. This is an admission check, not a proof that every accepted
calculation will fit: basis size, method, engine implementation, and other
running programs still affect peak memory.

## Results folder

Choose **Results folder** from File or Execution to place projects in a visible,
writable directory. The default on macOS and Linux is:

```text
~/Documents/OQP Studio/jobs
```

If that location cannot be written, Studio tries `~/OQP Studio/jobs`, its
application-support directory, and finally a temporary directory. The active
location is displayed in Execution. Choosing an explicit durable folder is
recommended because a system temporary directory may be cleaned automatically.

Each job directory contains its input, calculation log, metadata, and any JSON,
Molden, trajectory, cube, or other files written by OpenQP. These are ordinary
user files and can be copied to another computer or opened independently.

## Run and monitor

Select **Run** and follow **Live log**. This panel follows the actual calculation
log rather than the backend process console. Job states include queued,
running, cancelling, cancelled, done, not converged, and failed.

A result can be scientifically incomplete even when the process exits normally.
Studio marks a recognized unconverged optimization as **not converged**, retains
the best geometry, and offers **Restart** in Analysis. Restart prepares a new
input using that retained geometry and preserves the original runner and thread
selection when they remain available.

Use **Stop** on a queued or running job to request termination. Studio first
asks the process to terminate and escalates if it does not stop. Partial output
already written remains in the project directory for inspection.

## Existing calculations

Analysis can import files produced outside Studio. Select **Open existing
results**, then choose related log, JSON, Molden, trajectory, cube, input, and
coordinate files together. Studio copies them into a managed project so its
analysis tools can refer to one stable set of files.
