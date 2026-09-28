# examples/dnp_sol/solid_effect_freq_scan_2.m

- Signature: `solid_effect_freq_scan_2()`

## Purpose

Computes a powder-averaged, laboratory-frame DNP steady state for a 15N-labelled urea molecule at a specified distance from one electron, scanning two microwave-frequency windows. The source estimates a calculation time of hours.

## Model and method

The spin system contains one electron, two 15N nuclei, and four 1H nuclei at the coordinates defined in the source, in a 3.4 T field. It uses the `sphten-liouv` formalism, `IK-0` approximation, four-spin-order inter-level restriction, and projections [-2, -1, 0, +1, +2]. The secular Weizmann relaxation model has zero equilibrium state and temperature 4.2; electron/nuclear rates and the 7-by-7 distance-dependent R1d and R2d matrices are explicitly assigned in the source.

The electron is driven at 100 kHz. The scan concatenates 100 points from 144.0–145.5 MHz with 100 points from 14.0–15.5 MHz. Powder averaging uses the `rep_2ang_100pts_sph` grid, the `lvn-backs` method, and the ESR context in a call to `powder` with `dnp_freq_scan`.

## Outputs

The source plots the real longitudinal 1H and 15N expectation values against each frequency window, using the two columns of the result for the respective nuclear observables.
