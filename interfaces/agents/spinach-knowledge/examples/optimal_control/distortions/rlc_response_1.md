# examples/optimal_control/distortions/rlc_response_1.m

- Signature: `rlc_response_1(interp_type)`

## Purpose

Illustrates how a resonator response affects a typical composite pulse in NMR. `interp_type` defaults to `'previous'`, corresponding to a piecewise-constant input waveform; other `interp1` interpolation options, including `'linear'` and `'cubic'`, may also be used. Calculation time: seconds.

## Method

For ¹⁴N NMR at 14.09 T, the script creates random amplitude and phase controls over 20 slices spanning 50 μs, interpolates them onto a time grid set at four times the Nyquist rate, and forms the modulated input signal. It applies a second-order RLC band-pass response with Q = 80, heterodynes the output, and plots the input and output in both wall-clock and rotating-frame representations.
