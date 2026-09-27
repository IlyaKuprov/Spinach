# experiments/nmr_liquids/ecosy.m

- Signature: `fid=ecosy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive E.COSY pulse sequence. The implementation includes coherence orders through six-quantum order.

## Implementation

The kernel forms the initial `Lx` state and detection state, then propagates the F1 trajectory. It creates separate States-quadrature branches with `Lx` and `Ly` second pulses. The multiple-quantum filter combines coherence-order projections at orders `+/-2` through `+/-6`, with respective weights `1`, `2`, `4`, `6`, and `9`. A third pulse is applied to each branch and F2 detection returns real and imaginary components in `fid.cos` and `fid.sin`. The Liouvillian is `H+1i*R+1i*K` and the dwell time is `1/parameters.sweep`.

## Parameters / inputs

- `parameters.sweep` — sweep width in Hz.
- `parameters.npoints` — number of points for both dimensions.
- `parameters.spins` — nuclei on which the sequence runs, specified as `{'1H'}`, `{'13C'}`, etc.
- `H` — Hamiltonian matrix, received from context function.
- `R` — relaxation superoperator, received from context function.
- `K` — kinetics superoperator, received from context function.

The implementation requires `sphten-liouv` formalism, same-sized matrix inputs `H`, `R`, and `K`, a positive scalar sweep width, two positive integer point counts, and one isotope present in the system.

## Outputs

- `fid.cos`, `fid.sin` — real and imaginary components of the States quadrature signal.

## References

- [E.COSY reference](https://doi.org/10.1021/ja00308a042)
- [E.COSY reference](https://doi.org/10.1063/1.451421)
- [E.COSY reference](https://doi.org/10.1016/0022-2364(87)90102-8)
- [Spin Dynamics Wiki: `ecosy.m`](https://spindynamics.org/wiki/index.php?title=ecosy.m)
