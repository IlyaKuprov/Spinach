# experiments/nmr_liquids/pansy_triple.m

- Signature: `fid=pansy_triple(spin_system,parameters,H,R,K)`

## Purpose

Triple-channel PANSY pulse sequence. The source cites [10.1021/ja0634876](https://doi.org/10.1021/ja0634876) and [10.1007/128_2011_226](https://doi.org/10.1007/128_2011_226).

## Sequence and output

The routine forms `L = H + 1i*R + 1i*K` and starts from `Lz` magnetization on the first working spin. Positive- and negative-phase branches are created, +1 coherence is selected, and both branches evolve during F1. A second simultaneous pulse on all three working spins is followed by branch subtraction for axial-peak elimination. Detection on each spin yields:

- `fid.aa`: COSY FID on F1 nuclei.
- `fid.ab`: COSY FID on F1,F2 nuclei.
- `fid.ac`: COSY FID on F1,F3 nuclei.

Decoupling any working nucleus is not supported. This is the magnitude-mode analytical-pathway version; gradient echo/anti-echo PANSY variants are separate pulse-sequence functions. The routine requires `sphten-liouv` formalism.

## Inputs

- `parameters.spins`: three working nuclei, e.g. `{'1H','13C','15N'}`.
- `parameters.sweep`: three positive sweep widths in Hz.
- `parameters.npoints`: three positive integer point counts.
- `H`, `R`, and `K`: Hamiltonian matrix, relaxation superoperator, and kinetics superoperator from the context function; the matrices must have matching dimensions.

## Reference link

[Spinach Wiki: pansy_triple.m](https://spindynamics.org/wiki/index.php?title=pansy_triple.m)
