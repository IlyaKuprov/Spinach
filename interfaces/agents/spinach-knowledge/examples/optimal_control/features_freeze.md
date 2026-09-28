# examples/optimal_control/features_freeze.m

- Signature: `features_freeze()`

## Purpose

Demonstrates pulse optimisation with selected waveform samples held fixed while the remaining samples are adjusted. The example studies state-to-state transfer in a scalar-coupled hydrofluorocarbon spin system, transferring Z magnetisation from ¹H to ¹⁹F, and evaluates the final transfer fidelity.

## Physical / mathematical content

The control problem includes offset and RF-power sampling. Two pulse-sample ranges, 30–40 and 70–80, are frozen during optimisation; the other samples remain available to the optimiser.

## Numerical / algorithmic content

The pulse is optimised with Newton–Raphson GRAPE while enforcing the frozen-sample constraints. Starting from a random pulse, the optimisation typically reaches a fidelity of 0.999999, indicating near-unity transfer. The script then simulates the resulting waveform and reports its final fidelity. The source sets up six control channels.

## Implementation structure

The MATLAB code constructs the spin system, initial and target states, and controls; defines the frozen ranges; performs pulse optimisation over the offset and power sampling; and runs a final fidelity simulation. Source reference: http://dx.doi.org/10.1063/1.4949534
