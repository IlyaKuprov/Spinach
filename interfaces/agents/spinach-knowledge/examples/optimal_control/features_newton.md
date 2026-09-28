# examples/optimal_control/features_newton.m

- Signature: `features_newton()`

## Purpose

Optimise a pulse for state-to-state transfer across two scalar couplings in a three-spin hydrofluorocarbon fragment, from Lz on 1H to Lz on 19F. The six-channel pulse is designed to be robust to proton transmitter offset and reduced pulse nutation frequency. It uses Newton-Raphson GRAPE as described in [Goodwin and Kuprov (2016)](http://dx.doi.org/10.1063/1.4949534), with point-by-point waveform variation and a penalty when the waveform exceeds a user-specified power threshold. The initial guess is random. Calculation time: minutes.

## Physical / mathematical content

- The spin system contains 1H, 13C and 19F at a magnetic field of 9.4 T. All chemical shifts are 0 ppm; the 1H–13C and 13C–19F scalar couplings are 140 Hz and −160 Hz, respectively.
- The initial and target states are normalised Lz states on 1H and 19F. The six control operators are Lx and Ly for each isotope; the proton Lz operator defines the offset variation.
- The optimisation samples five proton offsets from −1 to +1 kHz and three pulse-power levels, 2π × [0.8, 0.9, 1.0] × 10³ rad/s, to account for nutation-frequency variation.

## Numerical / algorithmic content

- Uses the `sphten-liouv` formalism with no basis approximation. The drift Hamiltonian is constructed under the `nmr` assumption.
- Configures `control.method='newton'`, the `NS` penalty with weight 0.01, and a maximum of 50 iterations. The waveform has 100 slices of 10⁻⁴ s each; its initial guess is `randn(6,100)/10`.
- Calls `fmaxnewton(spin_system,@grape_xy,guess)`, scales the resulting pulse by the mean power level, and tests it with `shaped_pulse_xy` using `expv-pwc` propagation. The reported fidelity is `real(rho_targ'*rho)`.

## Implementation structure

1. Define the field, isotopes, chemical shifts, scalar couplings and basis; create the Spinach spin system.
2. Construct and normalise the initial and target states, then assemble the six control operators, proton offset operator and drift Hamiltonian.
3. Configure the control channels, robustness samples, slice durations, penalty, Newton optimisation and diagnostic plots; initialise the optimisation context with `optimcon`.
4. Generate a random initial waveform, optimise it, and simulate the resulting pulse to report the final target-state fidelity.
