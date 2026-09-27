# experiments/nmr_liquids/roesy.m

- Signature: `fid=roesy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive homonuclear ROESY with an ideal spin-lock. The source cites [10.1021/ja00315a069](https://doi.org/10.1021/ja00315a069) and [10.1016/0022-2364(85)90171-4](https://doi.org/10.1016/0022-2364(85)90171-4).

## Sequence and output

The routine forms `L = H + 1i*R + 1i*K`, applies a 90° pulse, and evolves during F1. The analytical spin-lock creates cosine and sine branches; each mixes under `1i*R + 1i*K` for `parameters.tmix`, then evolves and is detected during F2.

- Output: `fid.cos` and `fid.sin`, the free-induction decay components for hypercomplex processing.
- The ideal spin-lock model does not represent finite RF amplitude, RF offset, Hartmann-Hahn matching errors, or explicit RF phase transients.
- The routine accepts `sphten-liouv` and `zeeman-liouv` formalisms.

## Inputs

- `parameters.sweep`: two positive sweep widths in Hz, ordered F1,F2.
- `parameters.npoints`: two positive integer point counts, ordered F1,F2.
- `parameters.spins`: one working-spin label, e.g. `{'1H'}`.
- `parameters.tmix`: non-negative mixing time in seconds.
- `parameters.rho0`: initial state.
- `H`, `R`, and `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator from the context function; the matrices must have matching dimensions.

## Reference link

[Spinach Wiki: roesy.m](https://spindynamics.org/wiki/index.php?title=roesy.m)
