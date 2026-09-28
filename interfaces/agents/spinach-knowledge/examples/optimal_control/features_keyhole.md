# examples/optimal_control/features_keyhole.m

- Signature: `features_keyhole()`

## Purpose

Demonstrates a keyhole objective in optimal-control pulse design: selected two-spin-correlation terms are specified at an intermediate point in the pulse sequence, while the full pulse is optimised for the target state.

## Physical / mathematical content

The spin system provides 1H, 13C, and 19F channels with x- and y-phase controls. The keyhole constrains two-spin correlations at interval 20 in a 50-interval pulse.

## Numerical / algorithmic content

The objective is optimised with LBFGS-GRAPE over five RF-power levels. The numerical setup uses six controls (x and y for each of the three nuclei) and a 50-interval waveform.

## Implementation structure

The script builds and normalises the initial and target states, configures the keyhole correlation constraint and RF-power ensemble, optimises the waveform, and simulates the shaped pulse. It reports the target-state overlap after the final simulation.
