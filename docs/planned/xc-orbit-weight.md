# Plumb the existing XC atom-orbit weight into the gradient and response kernels

## The finding

`symAtomWeight` — the symmetry-orbit weight for the DFT quadrature — exists and
is consumed at `source/dftlib/dft_gridint.F90:2467` and `:2665`, but it has
**exactly one setter in the whole tree**: `dft_gridint_energy.F90:358`.

It is absent from:
- `source/dftlib/dft_gridint_grad.F90`
- `source/dftlib/dft_gridint_fxc.F90`
- `source/dftlib/dft_gridint_gxc.F90`
- `source/dftlib/dft_gridint_giao.F90`

So the reduction is written and validated on the energy path, and every gradient
recomputes the full grid.

## Value

1.5–3.0x on the **XC portion** of every DFT gradient — which means every
optimisation step, every NAMD step, every excited-state gradient. This is the
cheapest of the XC items because only the plumbing is missing.

Be skeptical about the response-kernel version: prior profiling on this machine
puts int2 METC at 96–98% of CPU MRSF response cost, so the `fxc`/`gxc` variants
are worth much less than they look.

## Verification

`E_xc` and `N_elec` against a C1 run to 1e-10; gradient parity to 1e-9.

## Machinery

Partial — the weight array exists; the consumers need the optional argument.
