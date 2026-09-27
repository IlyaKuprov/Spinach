# experiments/nmr_liquids/mqs_refocus.m

- Signature: `fid=mqs_refocus(spin_system,parameters,H,R,K)`

## Purpose

Two-dimensional multiple-quantum correlation with refocusing. The source cites [10.1063/1.432450](https://doi.org/10.1063/1.432450) and [10.1016/0022-2364(80)90096-7](https://doi.org/10.1016/0022-2364(80)90096-7).

## Sequence and output

The routine forms `L = H + 1i*R + 1i*K`, applies an initial 90° pulse, evolves for `parameters.delay_1`, applies a 180° refocusing pulse, then evolves for `parameters.delay_1` again. It selects the first requested coherence order, records the F1 trajectory, applies the final pulse of angle `parameters.angle`, and selects the second requested order. A `parameters.delay_2` period is followed by a 180° refocusing pulse and another `parameters.delay_2` period before F2 acquisition.

- Output: `fid`, a two-dimensional free-induction decay for amplitude-mode processing.
- In practice this implementation is homonuclear and uses exact analytical coherence-order projection rather than an explicit phase cycle.
- It requires the `sphten-liouv` formalism.

## Inputs

- `parameters.sweep`: two sweep widths in Hz, ordered [F1 F2].
- `parameters.npoints`: two point counts, ordered [F1 F2].
- `parameters.spins`: two spin labels; the documented example is `{'1H','1H'}`.
- `parameters.angle`: final-pulse flip angle in radians.
- `parameters.mqorder`: two integer coherence orders.
- `parameters.delay_1`, `parameters.delay_2`: evolution delays in seconds.
- `parameters.rho0`: initial state; `parameters.coil`: detection state.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics superoperators supplied by the context function; they must be matrices of matching dimensions.

## Reference link

[Spinach Wiki: mqs_refocus.m](https://spindynamics.org/wiki/index.php?title=mqs_refocus.m)
