# examples/fundamentals/convention_tests/nqi_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/convention_tests/nqi_test.m`
- Signature: `nqi_test()`
- Total lines: 39

## Purpose

Test of the reverse decomposition of spin-1 Hamiltonians.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Test of the reverse decomposition of
- spin-1 Hamiltonians.
- Get random test Hamiltonian
- Translate back
- Set up Spinach
- Re-build using Spinach functionality
- Compare the matrices

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `remtrace()`, `ham2nqi()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `orientation()`, `omega()`, `operator()`.
