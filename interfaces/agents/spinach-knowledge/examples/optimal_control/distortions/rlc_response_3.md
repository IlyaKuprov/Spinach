# examples/optimal_control/distortions/rlc_response_3.m

- Signature: `rlc_response_3()`

## Purpose

Probes how circuit response affects a deuterium pre-phasing pulse for the CD₃ group of alanine, intended to set the deuterium magnetisation for rephasing 100 μs after the pulse. The ensemble is a powder of 100 orientations with a B₁ distribution from 40 to 60 kHz per channel. It uses a piecewise-linear GRAPE pulse. Calculation time: minutes.

## Method

The deuterium pulse uses LBFGS with at most 100 iterations and 75 slices of 2 μs each. Its piecewise-linear waveform is passed through the resonator-response calculation for a probe circuit with Q = 200; the example assumes a 600 MHz magnet.
