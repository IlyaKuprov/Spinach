# examples/optimal_control/state_transfer_wf.m

- Signature: `state_transfer_wf()`

## Purpose

Transfer population from the lowermost to the uppermost energy level of a four-spin system using wave-function-space GRAPE. Calculation time: minutes.

## Physical / mathematical content

- The system contains two `1H` and two `13C` spins at a magnetic field of 14.1, with scalar Zeeman values `[1.5, 2.0, 30.0, 40.0]` and specified pairwise scalar couplings.
- The calculation uses the `zeeman-wavef` formalism without a basis approximation. The initial and target states are, respectively, the last and first entries of 16-element state vectors.
- The drift Hamiltonian includes transmitter offsets of 1050 for `1H` and 5285 for `13C`; the four controls are the `Lx` and `Ly` operators for each isotope.

## Numerical / algorithmic content

- The controls use 100 slices of duration `1.5e-4`, with pulse-power levels `2*pi*[460 480 500 520 540]`. The `NS` and `SNS` penalties have weights `[0.1 10]`.
- Starting from a random 4-by-100 waveform, `fmaxnewton` optimises `@grape_xy` with method `goodwin` and a maximum of 100 iterations. The resulting pulse is scaled by the mean power level.
- A test simulation applies the pulse with `shaped_pulse_xy` using `expv-pwc` propagation and reports `real(rho_targ'*rho)`.

## Implementation structure

- Create the spin system and basis, then construct the initial and target states, control operators, and offset-adjusted drift Hamiltonian.
- Configure the controls and optimisation, optimise the waveform, scale the pulse, and run the test simulation.
