# tests/kernel/test_dynamic_state_equilibrium_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_state_equilibrium_suite.m`
- Signature: `result=test_dynamic_state_equilibrium_suite()`
- Total lines: 73

## Purpose

Tests thermal equilibrium state construction paths. Syntax: result=test_dynamic_state_equilibrium_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Thermal equilibrium constructors\n')`.
- Lines 19-22: State the equilibrium-constructor target of the test; implemented by `result=new_test_result('kernel/dynamic_state_equilibrium_suite', 'Thermal equilibrium constructors', 'equilibrium() must produce normalised Boltzmann states in supported…`.
- Lines 24-25: Build a one-spin Hilbert-space system with finite temperature; implemented by `sys.magnet=0`.
- Lines 33-34: Set an explicit non-degenerate Hamiltonian in angular frequency units; implemented by `H=2*pi*diag([-5e11 3e11])`.
- Lines 39-40: Compare Hilbert-space equilibrium to the direct Boltzmann density matrix; implemented by `rho_h=equilibrium(spin_h,H)`.
- Lines 48-49: Build the matching Zeeman-Liouville system; implemented by `bas.formalism='zeeman-liouv'`.
- Lines 52-53: Compare left-product Liouville equilibrium to the vectorised density matrix; implemented by `H_left=kron(speye(2),H)`.
- Lines 58-59: Build a minimal anisotropic Hamiltonian cell for the orientation branch; implemented by `Q{1}=cell(3,3)`.
- Lines 65-66: Compare the four-argument branch to the explicitly oriented Hamiltonian; implemented by `rho_orient=equilibrium(spin_h,H,Q,euler_angles)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_state_equilibrium_suite', 'Thermal equilibrium constructors', 'equilibrium() must produce normalised Boltzmann states in supported…`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `inter.temperature` using `inter.temperature=300`.
- Lines 29: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_h` using `spin_h=test_spin_system(sys,inter,bas)`.
- Lines 34: computes `H` using `H=2*pi*diag([-5e11 3e11])`.
- Lines 35: computes `beta` using `beta=spin_h.tols.hbar/(spin_h.tols.kbol*spin_h.rlx.temperature)`.
- Lines 36: computes `rho_ref` using `rho_ref=expm(-beta*H)`.
- Lines 40: computes `rho_h` using `rho_h=equilibrium(spin_h,H)`.
- Lines 50: computes `spin_l` using `spin_l=test_spin_system(sys,inter,bas)`.
- Lines 53: computes `H_left` using `H_left=kron(speye(2),H)`.
- Lines 54: computes `rho_l` using `rho_l=equilibrium(spin_l,H_left)`.
- Lines 59: computes `Q{1}` using `Q{1}=cell(3,3)`.
- Lines 60: computes `[Q{1}{:}]` using `[Q{1}{:}]=deal(sparse(2,2))`.
- Lines 61: computes `Q{1}{2,2}` using `Q{1}{2,2}=sparse(2*pi*diag([2e11 -2e11]))`.

## Outputs

- result -regression test result with explanatory messages
- The test checks Hilbert-space, Zeeman-Liouville, and oriented-Hamiltonian
- thermal equilibrium construction against direct Boltzmann references.

## Implementation structure

- Tests thermal equilibrium state construction paths. Syntax:
- result=test_dynamic_state_equilibrium_suite()
- result -regression test result with explanatory messages
- The test checks Hilbert-space, Zeeman-Liouville, and oriented-Hamiltonian
- thermal equilibrium construction against direct Boltzmann references.
- Announce the test target
- State the equilibrium-constructor target of the test
- Build a one-spin Hilbert-space system with finite temperature
- Set an explicit non-degenerate Hamiltonian in angular frequency units
- Compare Hilbert-space equilibrium to the direct Boltzmann density matrix
- Build the matching Zeeman-Liouville system
- Compare left-product Liouville equilibrium to the vectorised density matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `equilibrium()`, `test_spin_system()`, `test_close()`, `speye()`, `rho_ref()`, `deal()`, `orientation()`.
