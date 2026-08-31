# tests/kernel/test_dynamic_remaining_spectral_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_remaining_spectral_suite.m`
- Signature: `result=test_dynamic_remaining_spectral_suite()`
- Total lines: 226

## Purpose

Tests remaining spectral, symmetry, and Fokker-Planck utilities. Syntax: result=test_dynamic_remaining_spectral_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file also defines local helper function(s): `local_created_system()`, `local_minimal_system()`, `local_tiny_rank_one()`, `gcp()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Remaining spectral and spatial utilities\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/dynamic_remaining_spectral_suite', 'Remaining spectral and spatial utilities', 'Small spectral, symmetry, and Fokker-Planck helpers must m…`.
- Lines 25-26: Check exact diagonalisation path in rspt_eig against eig(); implemented by `spin_system=local_minimal_system('zeeman-hilb',2)`.
- Lines 49-50: Check Liouville-space resonance-field extraction for a diagonal pencil; implemented by `spin_system=local_minimal_system('sphten-liouv',2)`.
- Lines 78-79: Check two-root Hilbert-space resonance extraction in a curved level gap; implemented by `spin_system=local_minimal_system('zeeman-hilb',2)`.
- Lines 106-107: Check one-dimensional gradient operator construction in Fokker-Planck space; implemented by `spin_system=local_created_system('zeeman-hilb',1)`.
- Lines 118-119: Check one-dimensional velocity generator against the finite-difference derivative; implemented by `spin_system=local_minimal_system('sphten-liouv',1)`.
- Lines 132-133: Check rotor-stack rank-zero assembly against direct Hamiltonian generation; implemented by `spin_system=local_created_system('sphten-liouv',14.1)`.
- Lines 151-152: Check S2 fully symmetric projector construction on a four-state product basis; implemented by `spin_system=local_minimal_system('sphten-liouv',4)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_remaining_spectral_suite', 'Remaining spectral and spatial utilities', 'Small spectral, symmetry, and Fokker-Planck helpers must m…`.
- Lines 26: computes `spin_system` using `spin_system=local_minimal_system('zeeman-hilb',2)`.
- Lines 27: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.
- Lines 28: computes `parameters.rho0` using `parameters.rho0=diag([0.9 0.1])`.
- Lines 29: computes `Hz` using `Hz=diag([1 3])`.
- Lines 30: computes `Hc` using `Hc=[0 0.1;0.1 0]`.
- Lines 31: computes `Hmw` using `Hmw=[0 1;1 0]`.
- Lines 32: computes `field` using `field=2`.
- Lines 33: computes `[E_obs,~,dE_obs,T_obs,LP_obs]` using `[E_obs,~,dE_obs,T_obs,LP_obs]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,field)`.
- Lines 34: computes `[V_ref,E_ref]` using `[V_ref,E_ref]=eig(field*Hz+Hc,'vector')`.
- Lines 35: computes `[E_ref,sort_idx]` using `[E_ref,sort_idx]=sort(E_ref,'ascend')`.
- Lines 36: computes `V_ref` using `V_ref=V_ref(:,sort_idx)`.
- Lines 37: computes `dE_ref` using `dE_ref=real(diag(V_ref'*Hz*V_ref))`.
- Lines 38: computes `T_ref` using `T_ref=abs(V_ref'*Hmw*V_ref).^2`.
- Lines 39: computes `LP_ref` using `LP_ref=real(diag(V_ref'*parameters.rho0*V_ref))`.
- Lines 51: computes `parameters` using `parameters=struct()`.
- Lines 52: computes `parameters.mw_freq` using `parameters.mw_freq=100`.
- Lines 53: computes `parameters.window` using `parameters.window=[9 11]`.

### Local helper functions

- Line 175: `local_created_system()` — `function spin_system=local_created_system(formalism,magnet)`. Build a quiet one-proton Spinach object in the requested formalism
  - Representative operation: `sys.magnet=magnet`.
  - Representative operation: `sys.isotopes={'1H'}`.
- Line 191: `local_minimal_system()` — `function spin_system=local_minimal_system(formalism,dim)`. Create a quiet minimal descriptor for matrix-only helper calls
  - Representative operation: `spin_system.sys.output='hush'`.
  - Representative operation: `spin_system.sys.enable={}`.
- Line 208: `local_tiny_rank_one()` — `function Q=local_tiny_rank_one(dim)`. Create a rank-one rotational-basis cell that evaluates to a tiny zero surrogate
  - Representative operation: `Q=cell(1,1)`.
  - Representative operation: `Q{1}=cell(3,3)`.
- Line 218: `local_ensure_pool()` — `function local_ensure_pool()`. Start a one-worker process pool for compact parfor utilities
  - Representative operation: `current_pool=gcp('nocreate')`.
  - Representative operation: `if isempty(current_pool)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks field-swept eigensystem helpers, rotor-stack assembly,
- permutation-symmetry projectors, and one-dimensional Fokker-Planck
- operators against compact analytical references.

## Implementation structure

- Tests remaining spectral, symmetry, and Fokker-Planck utilities. Syntax:
- result=test_dynamic_remaining_spectral_suite()
- result -regression test result with explanatory messages
- The test checks field-swept eigensystem helpers, rotor-stack assembly,
- permutation-symmetry projectors, and one-dimensional Fokker-Planck
- operators against compact analytical references.
- Announce the test target
- State the utility target of the test
- Check exact diagonalisation path in rspt_eig against eig()
- Check Liouville-space resonance-field extraction for a diagonal pencil
- Check two-root Hilbert-space resonance extraction in a curved level gap
- Check one-dimensional gradient operator construction in Fokker-Planck space

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_minimal_system()`, `rspt_eig()`, `V_ref()`, `test_close()`, `local_tiny_rank_one()`, `orientation()`, `eigenfields()`, `spin()`, `test_true()`, `isequal()`, `local_created_system()`, `g2fplanck()`, `hamiltonian()`, `assume()`, `spdiags()`.
