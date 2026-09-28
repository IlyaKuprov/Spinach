# examples/optimal_control/features_dt_opt.m

- Signature: `features_dt_opt()`

## Purpose

Optimises the durations of six pulse slices subject to a fixed total duration, then compares the resulting frequency-domain spectra. The durations are constrained to be nonnegative. The starting composite pulse is `270(-x)360(x)90(y)270(-y)360(y)90(x)` from Fig. 3 of <https://doi.org/10.1016/0022-2364(83)90133-6>; the example seeks a slightly improved pulse at the same power and total duration.

## Physical / mathematical content

The waveform is represented by six consecutive pulse segments whose durations are the optimisation variables; the total pulse duration is held fixed. The example models 100 non-interacting spins equally spaced across the affected spectral range, 25 kHz either side. It propagates the pulse sequence, acquires signals, and compares their Fourier-transformed spectra.

## Numerical / algorithmic content

The objective gradient is computed with tgrape. MATLAB fmincon uses an L-BFGS approximation while adjusting the six constrained durations.

## Implementation structure

The script defines the pulse and duration constraints, evaluates the objective and tgrape gradients in fmincon, propagates the optimised pulses, acquires and Fourier-transforms the signals, and prints the resulting slice durations. The sum of the six nonnegative durations is fixed.
