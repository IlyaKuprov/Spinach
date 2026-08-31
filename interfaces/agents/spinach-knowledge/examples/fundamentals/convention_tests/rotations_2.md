# examples/fundamentals/convention_tests/rotations_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/convention_tests/rotations_2.m`
- Signature: `rotations_2()`
- Total lines: 80

## Purpose

A rotations test comparing the Hamiltonians for a manually rotated (at the interaction specification level) spin system with the Hamiltonian that has been rotated using Spinach operator rotation functionality.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Generate random matrices; implemented by `shift_tensor_a=randn(3,3)`.
- Lines 15-18: % Kernel level rotation; implemented by `sys.magnet=14.1`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: A pair of spins at a distance, A; implemented by `sys.isotopes={'1H','15N'}`.
- Lines 33-34: Spinach housekeeping, A; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Hamiltonian, A; implemented by `[H,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 41-44: % Input level rotation; implemented by `sys.magnet=14.1`.
- Lines 50-51: A pair of spins at a distance, B; implemented by `sys.isotopes={'1H','15N'}`.
- Lines 59-60: Spinach housekeeping, B; implemented by `spin_system=create(sys,inter)`.
- Lines 63-64: Hamiltonian, B; implemented by `[H,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 67-70: % Comparison; implemented by `resnorm=norm(H_A-H_B,1)`.
- Lines 69-70: Norm of the residual; implemented by `resnorm=norm(H_A-H_B,1)`.
- Lines 72-73: Display the diagnostics; implemented by `if resnorm>1e-3`.

### Control flow inferred from the code

- Line 73: conditional branch on `resnorm>1e-3`.

### Key state/data transformations

- Lines 11: computes `shift_tensor_a` using `shift_tensor_a=randn(3,3)`.
- Lines 12: computes `shift_tensor_b` using `shift_tensor_b=randn(3,3)`.
- Lines 13: computes `coupl_tensor_a` using `coupl_tensor_a=100*randn(3,3)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `sys.isotopes` using `sys.isotopes={'1H','15N'}`.
- Lines 26-27: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={shift_tensor_a shift_tensor_b}`.
- Lines 28: computes `inter.coordinates` using `inter.coordinates={[0.7 0.8 0.9]`.
- Lines 30: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 31: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=coupl_tensor_a`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `[H,Q]` using `[H,Q]=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 39: computes `H_A` using `H_A=H+orientation(Q,[1.0 2.0 3.0])`.
- Lines 52: computes `R` using `R=euler2dcm(1.0,2.0,3.0)`.
- Lines 65: computes `H_B` using `H_B=H+orientation(Q,[0.0 0.0 0.0])`.
- Lines 70: computes `resnorm` using `resnorm=norm(H_A-H_B,1)`.

## Implementation structure

- A rotations test comparing the Hamiltonians for a manually rotated (at
- the interaction specification level) spin system with the Hamiltonian
- that has been rotated using Spinach operator rotation functionality.
- Generate random matrices
- % Kernel level rotation
- Magnet field
- Basis set
- A pair of spins at a distance, A
- Spinach housekeeping, A
- Hamiltonian, A
- % Input level rotation
- A pair of spins at a distance, B

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `orientation()`, `euler2dcm()`, `num2str()`, `report()`.
