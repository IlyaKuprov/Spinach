# experiments/nmr_liquids/pansy_cosy.m

- Signature: `fid=pansy_cosy(spin_system,parameters,H,R,K)`

## Purpose

Magnitude-mode PANSY-COSY. The source cites [10.1021/ja0634876](https://doi.org/10.1021/ja0634876) and [10.1016/j.pnmrs.2021.03.001](https://doi.org/10.1016/j.pnmrs.2021.03.001).

## Sequence and output

The routine forms `L = H + 1i*R + 1i*K` and starts from `Lz` magnetization on the first working spin. It creates positive- and negative-phase first-pulse branches, selects +1 coherence, evolves in F1, then applies the second pulse to the first and second spins. Subtracting the branches removes axial peaks. Detection on the first and second spins gives:

- `fid.aa`: magnitude-mode COSY FID on F1,F1 nuclei.
- `fid.ab`: magnitude-mode COSY FID on F1,F2 nuclei.

Decoupling either working nucleus is not supported. This is the magnitude-mode analytical-pathway version; gradient echo/anti-echo PANSY-COSY variants are separate pulse-sequence functions. The routine requires `sphten-liouv` formalism.

## Inputs

- `parameters.spins`: two working nuclei, e.g. `{'1H','13C'}`.
- `parameters.sweep`: two positive sweep widths in Hz.
- `parameters.npoints`: two positive integer point counts.
- `H`, `R`, and `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator from the context function; the matrices must have matching dimensions.

## Reference link

[Spinach Wiki: pansy_cosy.m](https://spindynamics.org/wiki/index.php?title=pansy_cosy.m)
