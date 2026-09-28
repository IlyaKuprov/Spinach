# examples/fundamentals/state_tests/thermal_equilibrium_3.m

- Signature: `thermal_equilibrium_3()`

## Purpose

Test of the thermal equilibrium functionality against the textbook expressions for the Boltzmann populations.

## Physical / mathematical content

Uses an X-band field of 0.34 T and a trityl electron with two protons. The example specifies anisotropic Zeeman principal values and Euler orientations, Cartesian coordinates, and a spin temperature of 80 K.

## Numerical / algorithmic content

For each of `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`, it computes the equilibrium state and the three spins' `Lz` expectation values. These are checked against textbook values formed from `levelpop` Boltzmann populations; each relative difference must be at most `1e-3`.

## Implementation structure

The loop creates the basis and equilibrium state for each formalism, evaluates the three observables, and compares Spinach results with the population-based expressions.
