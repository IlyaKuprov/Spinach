# examples/fundamentals/state_tests/thermal_equilibrium_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/thermal_equilibrium_1.m`
- Signature: `thermal_equilibrium_1()`
- Total lines: 57

## Purpose

Observables at thermal equilibrium using the three formalisms supported by Spinach kernel, tested against known answers.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Observables at thermal equilibrium using the three formalisms
- supported by Spinach kernel, tested against known answers.
- Spin system parameters
- Preallocate the answers
- Get numerical equilibrium magnetisation
- Get analytical equilibrium magnetisations
- Display the answers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `eq_mags_spinach()`, `state()`, `levelpop()`, `eq_mags_textbook()`, `any()`.
