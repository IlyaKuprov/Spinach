# examples/fundamentals/state_tests/thermal_equilibrium_4.m

- Signature: `thermal_equilibrium_4()`

## Purpose

Test of the invariance of the thermal equilibrium state under the thermalised relaxation superoperator.

## Physical / mathematical content

The model is a four-19F spin system at 9.4 T with damp relaxation, a 40 K temperature, and a damping rate of 5.0. The source tests both `sphten-liouv` and `zeeman-liouv` with full (`labframe`) relaxation retention.

## Numerical / algorithmic content

For each formalism, the script computes the equilibrium state and relaxation superoperator, then thermalises the latter using both the Dibari-Levitt method (`dibari`) and the inhomogeneous master equation method (`IME`). In each case it checks `norm(Rt*rho_eq,2)` against `1e-9`.

## Implementation structure

Constructs the basis and equilibrium state, builds relaxation, applies each thermalisation method, and verifies that the resulting superoperator annihilates the equilibrium state to the stated tolerance.
