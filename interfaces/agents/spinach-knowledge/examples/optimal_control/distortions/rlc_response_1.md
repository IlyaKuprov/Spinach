# examples/optimal_control/distortions/rlc_response_1.m

- MATLAB implementation: [examples/optimal_control/distortions/rlc_response_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/rlc_response_1.m)

Source: [examples/optimal_control/distortions/rlc_response_1.m](../../../../../../examples/optimal_control/distortions/rlc_response_1.m)

- Signature: `rlc_response_1(interp_type)`

## Purpose

Illustrates how a second-order RLC band-pass model changes a synthetic NMR carrier waveform. The script is a response demonstration rather than an optimal-control calculation; no imported acquisition or hardware response is used.

## Synthetic pulse and settings

The example is set for ¹⁴N NMR at 14.09 T, with angular frequency set by `omega = 14.09*spin('14N')` and quality factor Q = 80. It draws 21 random amplitude controls and 21 random phase controls for 20 slices over 50 μs, setting the first two and last three values of each control array to zero. The controls are interpolated to a time grid with step `pi/(4*omega)`, described in the source as four times the Nyquist rate. The optional `interp_type` is passed to `interp1`; it defaults to `'previous'` (piecewise constant), and the source also gives `'linear'` and `'cubic'` as examples.

The real carrier input is amplitude times cos(omega t + phase). Its input rotating-frame components are amplitude times cos(phase) and amplitude times sin(phase). The band-pass transfer function encoded in `tf` is `H(s) = (s/(omega*Q)) / (s^2/omega^2 + s/(omega*Q) + 1)`; `lsim` applies it to the input over the same time grid.

## Observables

The four-panel figure plots the input carrier and simulated output in wall-clock form, plus the input and output components heterodyned at omega. The output quadratures are formed by multiplying by cos(omega t) and sin(omega t), applying `lowpass` with argument 0.1, and using factors +2 and −2, respectively. Time is displayed in microseconds and voltage in arbitrary units; the rotating-frame traces are not experimental measurements.

The source comment estimates a calculation time of seconds; no timing benchmark is performed here. Random controls are not seeded, so traces can vary between runs.
