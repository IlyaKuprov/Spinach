# examples/fundamentals/state_tests/thermal_equilibrium_5.m

- Signature: `thermal_equilibrium_5()`

## Purpose

Cross-formalism test of state recovery towards the thermodynamic equilibrium.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Uses full (`labframe`) relaxation retention in both Liouville formalisms. For this damp-only model, damping is added after retention, so the generator is unchanged from the former diagonal setting; Zeeman diagonal retention is not supported.

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
