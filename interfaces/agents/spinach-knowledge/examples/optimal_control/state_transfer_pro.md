# examples/optimal_control/state_transfer_pro.m

- Signature: `state_transfer_pro()`

## Purpose

Optimise a pulse for magnetisation transfer from H(N) to C(O) in a typical protein backbone spin system, using literature shifts and couplings. The optimisation spans pulse powers to emulate B1 inhomogeneity and offsets to account for imperfect transmitter placement. It uses LBFGS-GRAPE with point-by-point waveform variation and a pulse-amplitude penalty.

## Physical / mathematical content

- The spin system contains six nuclei at a 9.4 T magnetic field, with specified chemical shifts and scalar couplings.
- The initial and target states are normalised `Lz` states on H and C, respectively. Control operators act on the `1H`, `13C`, and `15N` channels.
- The optimisation uses three offsets from −100 to 100 Hz on each channel and five pulse-power levels from 0.9 to 1.1 kHz, converted to rad/s.

## Numerical / algorithmic content

- The basis uses the `sphten-liouv` formalism and `IK-0` approximation. The pulse has 500 slices of 40 µs each.
- `fmaxnewton` optimises the waveform with `grape_xy` using the `lbfgs` method, an `NS` penalty of weight 0.01, and a maximum of 500 iterations.
- A test simulation applies the optimised pulse with `shaped_pulse_xy` and reports the real overlap with the target state.

## Implementation structure

- Define the magnetic field, spin system, chemical shifts, scalar couplings, and basis; then create the Spinach spin system.
- Construct the initial and target states, control and offset operators, drift Hamiltonian, and transmitter offsets.
- Configure the control ensemble and initialise a waveform guess that starts with a `1H` pulse and ends with a `13C` pulse. Optimise, scale the result by the mean power level, and run the test simulation.
- Estimated calculation time: hours. Reference: http://dx.doi.org/10.1016/j.jmr.2011.07.023
