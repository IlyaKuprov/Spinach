# examples/fundamentals/state_tests/thermal_equilibrium_2.m

- Signature: `thermal_equilibrium_2()`

## Purpose

Thermal equilibrium states, using all the different formalisms supported by Spinach kernel.

## Physical / mathematical content

Computes finite-temperature equilibrium for an E8, 1H, 14N, 15N spin system at 14.1 T and 4.2 K. The source includes several very large scalar couplings, annotated 'Remove these couplings to get machine precision.'

## Numerical / algorithmic content

The equilibrium density operator is calculated with `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv` (no basis approximation). Liouville representations are converted to the Zeeman-Hilbert matrix form, with the spherical-tensor result also transformed and scaled by the product of spin multiplicities. The test compares the resulting matrices and fails if either 2-norm difference exceeds `1e-8`.

## Implementation structure

Builds the same spin system and equilibrium state in each formalism, performs the representation conversions, then checks cross-formalism agreement.
