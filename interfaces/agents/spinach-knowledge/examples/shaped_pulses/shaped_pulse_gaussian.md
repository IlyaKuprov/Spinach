# examples/shaped_pulses/shaped_pulse_gaussian.m

- Signature: `shaped_pulse_gaussian()`

## Purpose

Gaussian 90-degree pulse on a chain of 31 strongly coupled protons. Calculation time: seconds

## Physical / mathematical content

- At a 14.1 T magnetic field, the 31 `1H` spins have scalar Zeeman values spanning -4 to 4, with 10 Hz scalar couplings between neighboring spins.
- The Gaussian waveform is read from `gaussian_1000.pk`, calibrated for a 90-degree rotation, and applied with a 480 Hz frequency offset.

## Numerical / algorithmic content

- Uses the `sphten-liouv` formalism with the `IK-2` approximation, `scalar_couplings` connectivity, and proximity level 1.
- Samples the 0.015 s pulse at 80 points, converts its amplitudes and phases to X and Y controls, and propagates the state with `shaped_pulse_xy` using `expv-pwc`.
- Acquires 2048 points over a 5000 Hz sweep, applies exponential apodisation with parameter 6, then Fourier-transforms with 16384-point zero filling and plots the imaginary spectrum.

## Implementation structure

- Creates the spin system and basis, applies NMR assumptions, and builds the Hamiltonian and `1H` control and offset operators.
- Starts from a `1H` `Lz` state, calibrates and executes the Gaussian pulse, then acquires, processes, and plots the resulting signal.
