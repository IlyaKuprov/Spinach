# examples/optimal_control/magic_pulse_phase.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/magic_pulse_phase.m)

## Purpose

This phase-control example sets up a broadband 90° 13C excitation pulse for the simultaneous transfers Lz to Lx, Ly to Ly, and Lx to -Lz. The model contains 100 non-interacting 13C spins equally distributed from -100 to +100 ppm at 28.18 T. The source cites the magic-pulse paper at [doi:10.1016/j.jmr.2005.12.010](https://doi.org/10.1016/j.jmr.2005.12.010) and motivates tolerance to resonance offsets and RF power-calibration variation.

## Design objective and constraints

The source describes a 50 microsecond duration ceiling based on a worst-case 13C-1H coupling of about 200 Hz, but configures 60 one-microsecond phase intervals (60 microseconds total). Those are distinct source statements; the example does not explain how the configured grid relates to the stated ceiling. It samples ten RF nutation levels from 50 to 70 kHz while holding the amplitude profile fixed, and uses GRAPE phase control through `fmaxnewton` with `@grape_phase` and `control.method='lbfgs'`. The initial phase profile is random, scaled by pi/5. The configured iteration limit is 200.

The normalised Lx, Ly, and Lz states define the three transfer targets. As in the Cartesian counterpart, the basis is `sphten-liouv` with `IK-2` at proximity level 1; the source comment says that this retains complete single-spin bases and omits multi-spin orders in this case.

## Evaluation shown by the example

The optimised phase and fixed amplitude are converted to Cartesian x/y controls and applied to an initial Lz state using the piecewise-constant exponential propagator. The resulting 13C signal is acquired with an L+ coil, a 70 kHz sweep, 2,048 points, 16,384-point zero filling, and a ppm axis; Gaussian apodisation with parameter 10 precedes the shifted Fourier transform. The script plots the real spectrum and compares it with a conventional hard-pulse spectrum configured with `2*pi*60e3` power, 4.2 microsecond duration, phase pi/2, and rank 3.

The source defines the design and plotting, but does not supply an observed waveform or convergence result. Its calculation-time comment says minutes. Source contacts: ilya.kuprov@weizmann.ac.il and david.goodwin@inano.au.dk.
