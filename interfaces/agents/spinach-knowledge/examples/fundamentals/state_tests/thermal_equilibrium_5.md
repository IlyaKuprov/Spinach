# examples/fundamentals/state_tests/thermal_equilibrium_5.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/thermal_equilibrium_5.m`
- Signature: `thermal_equilibrium_5()`
- Total lines: 86

## Purpose

Cross-formalism test of state recovery towards the thermodynamic equilibrium.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Cross-formalism test of state recovery towards the
- thermodynamic equilibrium.
- Magnet field
- Isotopes
- Chemical shifts
- J-couplings
- Relaxation theory parameters
- Formalisms and methods to test
- Loop over formalisms
- Loop over methods
- Thermalisation method
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `operator()`, `step()`, `assume()`, `hamiltonian()`, `relaxation()`, `state()`, `evolution()`.
