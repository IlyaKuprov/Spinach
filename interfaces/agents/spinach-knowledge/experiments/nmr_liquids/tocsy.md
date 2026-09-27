# experiments/nmr_liquids/tocsy.m

- Signature: `fid=tocsy(spin_system,parameters,H,R,K)`

## Purpose

Amplitude-mode homonuclear TOCSY experiment using continuous spin-lock mixing. The homonuclear spin species is selected with `parameters.spins`, for example `{'1H'}` or `{'13C'}`. The sequence is described in the cited TOCSY papers; it is not an explicit MLEV, DIPSI, WALTZ, or clean-TOCSY composite pulse-train simulation.

## Physical / mathematical content

The sequence starts from `parameters.rho0`, applies a 90-degree x pulse, and evolves the two-dimensional indirect and detection periods under the full Liouvillian `L = H + iR + iK`. During the mixing time it propagates under x- and y-oriented spin-lock terms, `L + 2 pi lamp Lx` and `L + 2 pi lamp Ly`, producing cosine and sine components for States quadrature processing. Here `lamp` is the spin-lock power in Hz. Relaxation and kinetics therefore act during the spin lock as well as during the other evolution periods.

## Parameters / inputs

- `parameters.sweep`: two positive sweep widths in Hz, for F1 and F2.
- `parameters.npoints`: two positive integer point counts, for F1 and F2.
- `parameters.spins`: a one-element cell array naming the working isotope, e.g. `{'1H'}`.
- `parameters.tmix`: non-negative mixing time in seconds.
- `parameters.lamp`: positive spin-lock power in Hz.
- `parameters.rho0`: initial state.
- `H`: Hamiltonian matrix; `R`: relaxation superoperator; `K`: kinetics superoperator, all supplied by the context function and required to have matching dimensions.

## Outputs

Returns `fid.cos` and `fid.sin`, the cosine and sine signal components used for States quadrature processing.

## References

- [TOCSY paper](https://doi.org/10.1016/0022-2364(83)90226-3)
- [TOCSY paper](https://doi.org/10.1021/ja00295a052)
- [TOCSY paper](https://doi.org/10.1016/0022-2364(85)90018-6)
- [Spin Dynamics Wiki: tocsy.m](https://spindynamics.org/wiki/index.php?title=tocsy.m)
