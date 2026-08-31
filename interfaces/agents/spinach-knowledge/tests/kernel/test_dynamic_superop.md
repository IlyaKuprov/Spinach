# tests/kernel/test_dynamic_superop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_superop.m`
- Signature: `result=test_dynamic_superop()`
- Total lines: 74

## Purpose

Tests superop() sparse-XYZ spherical-tensor product operators. Syntax: result=test_dynamic_superop()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_xyz_to_sparse()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Spherical-tensor superoperator construction\n')`.
- Lines 19-22: State the superop target of the test; implemented by `result=new_test_result('kernel/dynamic_superop', 'Spherical-tensor superoperator construction', 'superop() must produce sparse XYZ product superoperators consistent with…`.
- Lines 24-25: Build a two-spin spherical-tensor Liouville-space system; implemented by `sys.magnet=0`.
- Lines 33-34: Check the inactive-spin shortcut returns the full identity; implemented by `A_unit=local_xyz_to_sparse(superop(spin_system,[0 0],'left'),matrix_dim)`.
- Lines 38-39: Build left, right, commutator, and anticommutator forms for Lz on spin 1; implemented by `A_left=local_xyz_to_sparse(superop(spin_system,[2 0],'left'),matrix_dim)`.
- Lines 44-46: Check algebraic side-product identities; implemented by `result=test_close(result,'superop comm identity',A_comm,A_left-A_right,1e-15,1e-15, 'a commutator superoperator must equal left multiplication minus right multiplication…`.
- Lines 50-51: Check the Lz commutator eigenvalues on irreducible tensor projections; implemented by `[~,m_proj]=lin2lm(spin_system.bas.basis(:,1))`.
- Lines 55-56: Exercise a two-spin product operator through the multi-spin path; implemented by `A_left=local_xyz_to_sparse(superop(spin_system,[2 2],'left'),matrix_dim)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_superop', 'Spherical-tensor superoperator construction', 'superop() must produce sparse XYZ product superoperators consistent with…`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0,0}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 31: computes `matrix_dim` using `matrix_dim=size(spin_system.bas.basis,1)`.
- Lines 34: computes `A_unit` using `A_unit=local_xyz_to_sparse(superop(spin_system,[0 0],'left'),matrix_dim)`.
- Lines 39: computes `A_left` using `A_left=local_xyz_to_sparse(superop(spin_system,[2 0],'left'),matrix_dim)`.
- Lines 40: computes `A_right` using `A_right=local_xyz_to_sparse(superop(spin_system,[2 0],'right'),matrix_dim)`.
- Lines 41: computes `A_comm` using `A_comm=local_xyz_to_sparse(superop(spin_system,[2 0],'comm'),matrix_dim)`.
- Lines 42: computes `A_acomm` using `A_acomm=local_xyz_to_sparse(superop(spin_system,[2 0],'acomm'),matrix_dim)`.
- Lines 51: computes `[~,m_proj]` using `[~,m_proj]=lin2lm(spin_system.bas.basis(:,1))`.

### Local helper functions

- Line 67: `local_xyz_to_sparse()` — `function A=local_xyz_to_sparse(xyz,matrix_dim)`. Convert Spinach XYZ sparse triples into Matlab sparse storage
  - Representative operation: `A=sparse(xyz(:,1),xyz(:,2),complex(xyz(:,3)),matrix_dim,matrix_dim)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks the unit-operator shortcut, direct commutator and
- anticommutator identities, and Lz spherical-tensor projection eigenvalues.

## Implementation structure

- Tests superop() sparse-XYZ spherical-tensor product operators. Syntax:
- result=test_dynamic_superop()
- result -regression test result with explanatory messages
- The test checks the unit-operator shortcut, direct commutator and
- anticommutator identities, and Lz spherical-tensor projection eigenvalues.
- Announce the test target
- State the superop target of the test
- Build a two-spin spherical-tensor Liouville-space system
- Check the inactive-spin shortcut returns the full identity
- Build left, right, commutator, and anticommutator forms for Lz on spin 1
- Check algebraic side-product identities
- Check the Lz commutator eigenvalues on irreducible tensor projections

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `superop()`, `test_spin_system()`, `local_xyz_to_sparse()`, `test_close()`, `speye()`, `lin2lm()`, `test_true()`, `nnz()`, `xyz()`, `complex()`.
