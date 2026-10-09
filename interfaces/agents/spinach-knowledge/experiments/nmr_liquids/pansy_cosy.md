# experiments/nmr_liquids/pansy_cosy.m

- Signature: `fid=pansy_cosy(spin_system,parameters,H,R,K)`

## Purpose and pathway

Magnitude-mode PANSY-COSY. The source begins with `Lz` magnetisation on the first working spin and forms `L = H + 1i*R + 1i*K`. It creates two first-pulse branches using `+pi/2` and `-pi/2` rotations about the first spin's `Lx`, selects the `+1` coherence on that spin, and evolves each branch during F1. A second `pi/2` pulse about the sum of the two working spins' `Ly` operators is applied; subtracting the branches performs the source's axial-peak elimination. No further coherence-order filter is applied in the function. The header notes that neither working nucleus can be decoupled.

## Inputs and acquisition

- `parameters.spins`: two working nuclei in a cell array; source example: `{'1H','13C'}`.
- `parameters.sweep`: two positive sweep widths in Hz; `parameters.npoints`: two positive integer point counts, ordered F1 then the two detected dimensions.
- `H`, `R`, and `K`: numeric, same-sized matrices supplied by the context function. The source requires the `sphten-liouv` formalism.

The F1 trajectory contains `npoints(1)` states. The output fields are two-dimensional FIDs: `fid.aa` detects with `L+` on spin 1 and has shape `npoints(1) × npoints(1)`; `fid.ab` detects with `L+` on spin 2 and has shape `npoints(2) × npoints(1)` (direct F2 samples in rows, F1 states in columns). These are the source's F1/F1 and F1/F2 magnitude-mode COSY channels.

## References

- [PANSY paper](https://doi.org/10.1021/ja0634876)
- [PANSY paper](https://doi.org/10.1016/j.pnmrs.2021.03.001)
- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/pansy_cosy.m)
- [Spinach Wiki: pansy_cosy.m](https://spindynamics.org/wiki/index.php?title=pansy_cosy.m)
