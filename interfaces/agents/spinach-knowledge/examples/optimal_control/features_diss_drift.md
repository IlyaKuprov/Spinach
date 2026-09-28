# examples/optimal_control/features_diss_drift.m

- Signature: `features_diss_drift()`

## Purpose

Demonstrates optimal-control pulse optimisation with a dissipative drift and robustness to offset and RF-power variation. The pulse is optimised with LBFGS-GRAPE and then evaluated in the specified model.

## Physical / mathematical content

The spin dynamics include a dissipative drift. The objective is evaluated over offset and power ensembles, and includes an RF-amplitude penalty.

## Numerical / algorithmic content

The source configures 500 pulse slices and uses the lbfgs method with GRAPE derivatives (fmaxnewton with grape_xy). This is a limited-memory quasi-Newton optimisation.

## Implementation structure

The script builds the spin system and control operators, specifies offset and power ensembles, optimises the 500-slice waveform with the amplitude penalty, and reports the resulting performance.
