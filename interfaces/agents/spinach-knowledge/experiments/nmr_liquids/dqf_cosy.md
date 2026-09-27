# experiments/nmr_liquids/dqf_cosy.m

- Signature: `fid=dqf_cosy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive double-quantum filtered COSY pulse sequence.

## Implementation

The sequence starts from `Lz` magnetisation on `parameters.spins{1}`, applies an `Lx` 90-degree pulse, and records the F1 trajectory. Separate `Lx` and `Ly` second-pulse branches form the States quadrature components. Each branch is filtered to coherence orders `+2` and `-2` by an exact analytical coherence-order projection, then receives a third `Lx` 90-degree pulse before F2 detection. The two observable evolutions are returned as `fid.cos` and `fid.sin`; the filter is not an explicit phase cycle or finite-gradient selection block. Evolution uses `L=H+1i*R+1i*K` and timestep `1/parameters.sweep`.

## Parameters / inputs

- `parameters.sweep` — sweep width in Hz.
- `parameters.npoints` — number of points for both dimensions.
- `parameters.spins` — nuclei on which the sequence runs, specified as `{'1H'}`, `{'13C'}`, etc.; the selected isotope must have at least two spins.
- `H` — Hamiltonian matrix, received from context function.
- `R` — relaxation superoperator, received from context function.
- `K` — kinetics superoperator, received from context function.

The implementation requires `sphten-liouv` formalism, same-sized matrix inputs `H`, `R`, and `K`, a positive scalar sweep width, and two positive integer point counts.

## Outputs

- `fid.cos`, `fid.sin` — components of the free induction decay for hypercomplex processing.

## References

- [Double-quantum-filtered COSY reference](https://doi.org/10.1016/0006-291X(83)91225-1)
- [COSY reference](https://doi.org/10.1021/ja00388a062)
- [Spin Dynamics Wiki: `dqf_cosy.m`](https://spindynamics.org/wiki/index.php?title=dqf_cosy.m)
