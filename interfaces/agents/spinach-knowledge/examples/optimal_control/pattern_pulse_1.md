# examples/optimal_control/pattern_pulse_1.m

- Signature: `pattern_pulse_1()`

## Purpose

Design a phase-modulated pulse for nutation-frequency-selective excitation, as described in the Glaser group paper (https://doi.org/10.1016/j.jmr.2004.12.005). The pulse drives magnetisation into specified target states over specified nutation-frequency intervals. Calculation time: minutes.

## Physical / mathematical content

- A single on-resonance `13C` spin is simulated at a magnetic field of 28.18 T, with no basis approximations. The drift Hamiltonian is zero.
- The initial state is normalised `Sz`. Across 128 nutation frequencies from 6 to 14 kHz, the target is normalised `Sz` in the first 20, middle 20 (indices 54–73), and last 20 frequency samples, and normalised `Sx` elsewhere.
- The optimised control is the phase of a pulse with a fixed amplitude profile. Its Cartesian components act through `Lx` and `Ly`.

## Numerical / algorithmic content

- The script runs GRAPE phase optimisation with the L-BFGS method (`control.method='lbfgs'`, `@grape_phase`), starting from a constant `pi/4` phase profile and allowing up to 200 iterations.
- The pulse has 250 intervals of 20 µs each. The ensemble uses B1 power levels of `2*pi*nutf_range` rad/s and a separate initial–target state pair for each level.
- After optimisation, a `parfor` loop simulates the shaped pulse at each power level using `shaped_pulse_xy` with the `expv-pwc` propagator. The resulting states are projected onto `Sx` and `Sz` and plotted against the target pattern.

## Implementation structure

- Initialise the spin system and obtain the normalised `Sx` and `Sz` states, `Lx` and `Ly` control operators, and drift Hamiltonian.
- Construct and plot the frequency-dependent `Sx` and `Sz` target pattern.
- Set the control ensemble, pulse grid, plotting options, and initial phase guess; configure optimisation with `optimcon` and run `fmaxnewton` using `@grape_phase`.
- Simulate the optimised pulse across the nutation-frequency range in parallel and plot its `Sx` and `Sz` projections alongside the targets.
