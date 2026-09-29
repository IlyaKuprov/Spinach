# examples/optimal_control/features_freeze.m

- Signature: `features_freeze()`
- Source: [examples/optimal_control/features_freeze.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_freeze.m)
- Source reference: http://dx.doi.org/10.1063/1.4949534

## Purpose and model

This example optimises a simulated 1H-to-19F longitudinal-magnetisation transfer in a scalar-coupled 1H-13C-19F model at 9.4 T. Chemical shifts are zero; the 1H-13C and 13C-19F couplings are 140 Hz and -160 Hz. No measured or imported data are used.

## Controls and frozen samples

Six x/y controls address the three nuclei. The waveform has 100 intervals of 100 microseconds each (10 ms total). Optimisation samples five 1H offsets from -1000 to +1000 Hz and three RF powers, `2*pi*[800 900 1000]` rad/s. Samples 30-40 and 70-80, inclusive, are marked frozen for every control. The initial guess is random except those entries are set to 0.1; the source comment annotates them as “Frozen at 100 Hz.”

The source configures the NS penalty (weight 0.01), `control.method='goodwin'`, and 100 iterations, then calls `fmaxnewton` with `@grape_xy` (described in the source comments as Newton-Raphson GRAPE). Correlation-order, per-spin, x/y-control, and spectrogram plots are enabled. The source comment reports a typical optimisation fidelity of 0.999999.

## Output and limits

The waveform is scaled by the mean power and propagated with `shaped_pulse_xy` using `expv-pwc`; the script prints the real target-state overlap. This final call evaluates one drift model and does not pass the sampled offset ensemble, so the printed overlap is not a per-offset or per-power report. The page documents a simulation example, not hardware validation.
