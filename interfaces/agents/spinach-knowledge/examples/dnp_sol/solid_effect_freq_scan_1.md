# examples/dnp_sol/solid_effect_freq_scan_1.m

- Signature: `solid_effect_freq_scan_1()`

## Purpose

Calculates a single-crystal, laboratory-frame DNP steady state for a 15N-labelled urea molecule at one specified orientation and distance from one electron, while scanning two microwave-frequency windows. The source describes the calculation as taking minutes.

## Model and method

The spin system contains one electron, two 15N nuclei, and four 1H nuclei. It uses a 3.4 T field and the `sphten-liouv` formalism with the `IK-0` approximation, four-spin-order inter-level restriction, and projection set [-2, -1, 0, +1, +2]. Relaxation is the secular Weizmann DNP model with zero equilibrium state, temperature 4.2, specified electron/nuclear longitudinal and transverse rates, and distance-dependent rates set to 10^-3 for both R1d and R2d matrices.

The electron is driven with a 100 kHz microwave amplitude. The scan concatenates 100 points from 144.0–145.5 MHz and 100 points from 14.0–15.5 MHz. A single orientation [pi/4, pi/5, pi/6] is passed to `crystal`; the calculation calls `dnp_freq_scan` with the ESR context.

## Outputs

Four plots show the real longitudinal expectation values for 1H and 15N across each of the two frequency windows. The source forms the coil observables from the 1H and 15N longitudinal operators and uses the corresponding first and second result columns for those plots.
