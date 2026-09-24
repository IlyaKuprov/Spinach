# examples/fundamentals/state_tests/thermal_equilibrium_4.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/thermal_equilibrium_4.m`
- Signature: `thermal_equilibrium_4()`
- Total lines: 72

## Purpose

Test of the invariance of the thermal equilibrium state under the thermalised relaxation superoperator.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Uses full (`labframe`) relaxation retention in both Liouville formalisms. For this damp-only model, damping is added after retention, so the generator is unchanged from the former diagonal setting; Zeeman diagonal retention is not supported.

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Test of the invariance of the thermal equilibrium state under the
- thermalised relaxation superoperator.
- Magnet field
- Isotopes
- Chemical shifts
- J-couplings
- Relaxation theory parameters
- Formalisms to test
- Loop over formalisms
- Basis set
- Spinach housekeeping
- Isotropic thermal equilibrium

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `relaxation()`, `hamiltonian()`, `assume()`, `thermalize()`.
