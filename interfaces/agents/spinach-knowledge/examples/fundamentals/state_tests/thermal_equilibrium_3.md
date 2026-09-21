# examples/fundamentals/state_tests/thermal_equilibrium_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/thermal_equilibrium_3.m`
- Signature: `thermal_equilibrium_3()`
- Total lines: 73

## Purpose

Test of the thermal equilibrium functionality against the textbook expressions for the Boltzmann populations.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Test of the thermal equilibrium functionality against the
- textbook expressions for the Boltzmann populations.
- X-band magnet
- Electron and two protons
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Cartesian coordinates
- Spin temperature
- Formalisms to test
- Loop over formalisms
- Formalism and basis set
- Spinach housekeeping
- Isotropic thermal equilibrium

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `state()`, `levelpop()`.
