# experiments/nmr_liquids/roesy.m

- Signature: `fid=roesy(spin_system,parameters,H,R,K)`

## Purpose and sequence

Phase-sensitive homonuclear ROESY under the source's ideal-spin-lock model. The routine forms `L = H + 1i*R + 1i*K`, applies a `pi/2` pulse about `Lx` on the working spin to `parameters.rho0`, and records its F1 trajectory. The analytical spin-lock operator generates separate cosine and sine branches. During `parameters.tmix`, each branch evolves under `1i*R + 1i*K` (not `H`); F2 acquisition then uses the full `L` and detects `L+` on the same spin. The code does not insert a coherence-order selection step. It returns `fid.cos` and `fid.sin` for hypercomplex processing, each with F2 time samples in rows and stacked F1 states in columns: `npoints(2) × npoints(1)`.

This idealised model does not represent finite RF amplitude, RF offset, Hartmann-Hahn matching errors, or explicit RF phase transients.

## Inputs

- `parameters.spins`: one working spin in a cell array; source example: `{'1H'}`.
- `parameters.sweep`: two positive sweep widths in Hz; `parameters.npoints`: two positive integer point counts, in F1/F2 order.
- `parameters.tmix`: non-negative mixing time in seconds; `parameters.rho0`: numeric initial state with a Liouville-space dimension matching `H`.
- `H`, `R`, and `K`: numeric, same-sized matrices from the context function. Supported formalisms are `sphten-liouv` and `zeeman-liouv`.

## References

- [ROESY paper](https://doi.org/10.1021/ja00315a069)
- [ROESY paper](https://doi.org/10.1016/0022-2364(85)90171-4)
- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/roesy.m)
- [Spinach Wiki: roesy.m](https://spindynamics.org/wiki/index.php?title=roesy.m)
