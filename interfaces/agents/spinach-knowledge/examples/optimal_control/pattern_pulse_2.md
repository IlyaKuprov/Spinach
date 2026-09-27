# examples/optimal_control/pattern_pulse_2.m

- Signature: `pattern_pulse_2()`

## Purpose

Transmitter offset selective excitation described in Glaser group paper (https://doi.org/10.1016/j.jmr.2004.12.005). User-specified transmitter offset intervals have magnetisation arriving into user-specified states. The pulse is phase-modulated. Calculation time: minutes.

## Physical / mathematical content

- The script optimises a phase-modulated pulse to drive a single `13C` spin from `Sz` towards offset-dependent `Sx` or `Sz` target states across 128 transmitter offsets from −4 to 4 kHz.
- The pulse has 300 intervals of 20 µs each, with a fixed amplitude profile and a power level of `2*pi*2000`; the phase profile is optimised using LBFGS GRAPE (`fmaxnewton` with `grape_phase`).

## Numerical / algorithmic content

- After optimisation, the script simulates the pulse at each transmitter offset in a `parfor` loop and plots the resulting `Sx` and `Sz` projections against the targets.

## Implementation structure

- Transmitter offset selective excitation described in Glaser group
- paper (https://doi.org/10.1016/j.jmr.2004.12.005). User-specified
- transmitter offset intervals have magnetisation arriving into user-specified states. The pulse is phase-modulated.
- Calculation time: minutes.
- Magnetic field
- Single carbon spin
- Put the spin at 0 ppm
- No approximations
- Run Spinach housekeeping
- Get pertinent spin states
- Get pertinent control operators
