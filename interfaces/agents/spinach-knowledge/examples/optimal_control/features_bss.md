# examples/optimal_control/features_bss.m

- Signature: `features_bss()`

## Purpose

Optimises a 90-degree pulse for a single proton with a 1 MHz Larmor frequency when Bloch–Siegert corrections are included. The control is a significant fraction of the Larmor frequency, so the counter-rotating field shifts the resonance. The example compares pulses optimised with the correction enabled and disabled, evaluating both in the corrected model. Its offset ensemble is deliberately applied through transverse Lx rather than the usual Lz operator.

## Physical / mathematical content

The single-spin pulse-design problem includes Bloch–Siegert shift corrections and an offset ensemble represented by a transverse Lx term. The target is a 90-degree rotation.

## Numerical / algorithmic content

The pulse is optimised with LBFGS-GRAPE, once with Bloch–Siegert corrections enabled and once without them. Both resulting pulses are evaluated using the corrected model, and their fidelities are reported.

## Implementation structure

The script creates the one-proton spin system, sets Lx and Ly as RF controls, defines the transverse offset ensemble, runs the two optimisation cases, and reports corrected-model fidelities. The stated Larmor frequency is 1 MHz.
