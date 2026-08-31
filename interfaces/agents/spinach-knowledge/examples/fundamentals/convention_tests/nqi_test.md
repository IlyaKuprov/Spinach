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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Get random test Hamiltonian; implemented by `H_T=randn(3,3)+1i*randn(3,3)`.
- Lines 12-13: Translate back; implemented by `[omega,Q]=ham2nqi(H_T)`.
- Lines 15-16: Set up Spinach; implemented by `sys.magnet=0`.
- Lines 24-25: Re-build using Spinach functionality; implemented by `spin_system=assume(spin_system,'labframe')`.
- Lines 32-33: Compare the matrices; implemented by `if norm(H_T-H_S,2)>1e-6*norm(H_T+H_S,2)`.

### Control flow inferred from the code

- Line 33: conditional branch on `norm(H_T-H_S,2)>1e-6*norm(H_T+H_S,2)`.

### Key state/data transformations

- Lines 9: computes `H_T` using `H_T=randn(3,3)+1i*randn(3,3)`.
- Lines 13: computes `[omega,Q]` using `[omega,Q]=ham2nqi(H_T)`.
- Lines 16: computes `sys.magnet` using `sys.magnet=0`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'14N'}`.
- Lines 18: computes `inter.coupling.matrix` using `inter.coupling.matrix={Q/(2*pi)}`.
- Lines 19: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 20: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 21: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 26: computes `[H_iso,H_aniso]` using `[H_iso,H_aniso]=hamiltonian(spin_system)`.
- Lines 27: computes `H_S` using `H_S=H_iso+orientation(H_aniso,[0 0 0])`.

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
