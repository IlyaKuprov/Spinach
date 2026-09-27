# experiments/nmr_liquids/noesy.m

- Signature: `fid=noesy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive homonuclear NOESY. The source cites [10.1063/1.438208](https://doi.org/10.1063/1.438208), [10.1016/0006-291X(80)90695-6](https://doi.org/10.1016/0006-291X(80)90695-6), and [10.1016/0022-2364(82)90279-7](https://doi.org/10.1016/0022-2364(82)90279-7).

## Sequence and output

The routine forms `L = H + 1i*R + 1i*K`, applies the first 90° pulse, evolves during F1, and runs a four-step phase cycle. By default, homospoil retains longitudinal magnetisation before mixing under relaxation and kinetics (`1i*R + 1i*K`). Setting `parameters.oldschool` true disables this homospoil path and uses the full generator during mixing. The four acquisitions are combined by axial-peak elimination.

- Output: `fid.cos` and `fid.sin`, the two FID components for hypercomplex F1 processing.
- The layout is optimized for memory rather than CPU time and is intended for very large protein and nucleic-acid simulations.
- Non-empty analytical decoupling is meaningful only in `sphten-liouv`; the routine also accepts `zeeman-liouv` formalism.

## Inputs

- `parameters.sweep`: two positive sweep widths in Hz.
- `parameters.npoints`: two positive integer point counts.
- `parameters.spins`: one working-spin label, e.g. `{'1H'}` or `{'13C'}`.
- `parameters.tmix`: non-negative mixing time in seconds.
- `parameters.decouple` (optional): spin labels such as `{'13C','1H'}` or a numeric list of spin indices.
- `parameters.rho0`: initial state; exact thermal equilibrium can be requested through `parameters.needs={'rho_eq'}`.
- `parameters.oldschool` (optional): logical scalar; true disables the default homospoil gradient before mixing.
- `H`, `R`, and `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator from the context function.

## Reference link

[Spinach Wiki: noesy.m](https://spindynamics.org/wiki/index.php?title=noesy.m)
