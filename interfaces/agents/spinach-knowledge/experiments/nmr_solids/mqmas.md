# experiments/nmr_solids/mqmas.m

- Signature: `fid=mqmas(spin_system,parameters,H,R,K)`

## Purpose

Rotor-synchronous MQMAS pulse sequence for a 2D amplitude-mode free induction decay. Call it from the `singlerot.m` context, which supplies `H`, `R`, and `K`. See <https://spindynamics.org/wiki/index.php?title=mqmas.m>.

## Parameters / inputs

- `spin_system`: Spinach spin system.
- `H`, `R`, `K`: Numeric square matrices of equal size, supplied by `singlerot.m`.
- `parameters.spins`: Cell array containing one isotope string present in the spin system; selects the active spin.
- `parameters.pulse_dur`: Durations of the two pulses, in seconds; a two-element vector of non-negative real numbers.
- `parameters.pulse_amp`: Amplitudes of the two pulses, in rad/s; a two-element vector of real numbers.
- `parameters.mq_order`: Integer MQMAS coherence order.
- `parameters.rho0`: Initial condition, usually `Lz`.
- `parameters.coil`: Detection state, usually `L+`.
- `parameters.spc_dim`: Positive integer spatial problem dimension.
- `parameters.npoints`: Two-element vector of positive integers specifying the point counts in the indirect and direct dimensions.
- `parameters.rate`: Non-zero real MAS rate.
- `parameters.sweep`: Positive real sweep width; must equal `abs(parameters.rate)`. Both dimensions are sampled stroboscopically at this sweep width, relative to the rotor period.
- `parameters.decouple`: Cell array of isotope strings to decouple, or an empty cell array. Listed isotopes must be present in the system; analytical decoupling requires the `sphten-liouv` formalism.
- Other parameters required by the `singlerot.m` context.

## Outputs

- `fid`: 2D amplitude-mode free induction decay.

## Implementation summary

The sequence forms `L=H+1i*R+1i*K`, applies decoupling, and runs the first pulse before selecting `parameters.mq_order` coherence. It evolves the indirect dimension at intervals of `1/abs(parameters.rate)`, runs the second pulse, selects `+1` coherence, and evolves the direct dimension while detecting with `parameters.coil`.