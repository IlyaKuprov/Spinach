# tests/kernel/test_optimcon_grape_one_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_optimcon_grape_one_spin.m`
- Signature: `result=test_optimcon_grape_one_spin()`
- Total lines: 147

## Purpose

Tests one-spin optimal-control setup and Hilbert-space GRAPE. Syntax: result=test_optimcon_grape_one_spin()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file also defines local helper function(s): `local_spin_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: One-spin optimal-control setup and GRAPE gradient\n')`.
- Lines 20-23: State the optimal-control target of the test; implemented by `result=new_test_result('optimcon/grape_one_spin', 'One-spin optimal-control setup and GRAPE gradient', 'optimcon() and grape_hilb() must handle a tiny Hilbert-space cont…`.
- Lines 25-26: Ensure that a parallel pool is available for the ensemble loop; implemented by `current_pool=gcp('nocreate')`.
- Lines 31-32: Build a minimal Hilbert-space Spinach object for setup and GRAPE; implemented by `spin_system=local_spin_system()`.
- Lines 34-35: Define one-spin operators and a two-step timing grid; implemented by `S=pauli(2)`.
- Lines 39-40: Configure a minimal optimal-control problem; implemented by `control.isotopes={'E'}`.
- Lines 57-59: Check that optimcon() absorbed the control metadata consistently; implemented by `result=test_true(result,'optimcon control count',spin_system.control.ncontrols==1, 'one supplied control operator must be registered')`.
- Lines 65-66: Check that distortions work when the optimiser method is defaulted; implemented by `control_default=rmfield(control,'method')`.
- Lines 75-76: Check that exact-Hessian methods remain incompatible with distortions; implemented by `control_hessian=control`.
- Lines 89-90: Evaluate the GRAPE fidelity and gradient for a non-trivial waveform; implemented by `waveform=[7 11]`.
- Lines 94-95: Build the independent exact Hilbert-space trajectory; implemented by `rho=S.x`.
- Lines 103-105: Check the fidelity against direct matrix exponentiation; implemented by `result=test_close(result,'grape_hilb fidelity',fidelity,fid_ref,1e-13,1e-13, 'GRAPE fidelity must match direct Hilbert-space propagation')`.
- Lines 107-108: Compute a centred finite-difference gradient of the GRAPE fidelity; implemented by `step_size=1e-6`.
- Lines 122-124: Check the analytic GRAPE gradient against finite differences; implemented by `result=test_close(result,'grape_hilb gradient',grad,grad_ref,1e-8,1e-8, 'GRAPE gradient must match centred finite differences')`.
- Lines 126-128: Check that no heavy trajectory is returned when plotting is disabled; implemented by `result=test_true(result,'grape_hilb trajectory suppression',isempty(traj_data.forward), 'empty plotting requests should suppress trajectory storage')`.

### Control flow inferred from the code

- Line 27: conditional branch on `isempty(current_pool)`.
- Line 96: `for` loop over `n=1:numel(pulse_dt)`.
- Line 110: `for` loop over `n=1:numel(waveform)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('optimcon/grape_one_spin', 'One-spin optimal-control setup and GRAPE gradient', 'optimcon() and grape_hilb() must handle a tiny Hilbert-space cont…`.
- Lines 26: computes `current_pool` using `current_pool=gcp('nocreate')`.
- Lines 32: computes `spin_system` using `spin_system=local_spin_system()`.
- Lines 35: computes `S` using `S=pauli(2)`.
- Lines 36: computes `drift` using `drift=sparse(2,2)`.
- Lines 37: computes `pulse_dt` using `pulse_dt=[0.02 0.03]`.
- Lines 40: computes `control.isotopes` using `control.isotopes={'E'}`.
- Lines 41: computes `control.channels` using `control.channels=1`.
- Lines 42: computes `control.operators` using `control.operators={S.y}`.
- Lines 43: computes `control.rho_init` using `control.rho_init={S.x}`.
- Lines 44: computes `control.rho_targ` using `control.rho_targ={S.z}`.
- Lines 45: computes `control.pwr_levels` using `control.pwr_levels=1`.
- Lines 46: computes `control.pulse_dt` using `control.pulse_dt=pulse_dt`.
- Lines 47: computes `control.drifts` using `control.drifts={{drift}}`.
- Lines 48: computes `control.method` using `control.method='lbfgs'`.
- Lines 49: computes `control.max_iter` using `control.max_iter=0`.
- Lines 50: computes `control.penalties` using `control.penalties={'none'}`.
- Lines 51: computes `control.p_weights` using `control.p_weights=0`.

### Local helper functions

- Line 133: `local_spin_system()` — `function spin_system=local_spin_system()`. Build a minimal quiet Spinach object for Hilbert-space helpers
  - Representative operation: `spin_system.sys.output='hush'`.
  - Representative operation: `spin_system.sys.enable={}`.

## Outputs

- result -regression test result with explanatory messages
- The test checks that optimcon() accepts a minimal one-spin Hilbert-space
- control problem, and that grape_hilb() fidelity and gradient agree with
- independent matrix exponentiation and finite differences.

## Implementation structure

- Tests one-spin optimal-control setup and Hilbert-space GRAPE. Syntax:
- result=test_optimcon_grape_one_spin()
- result -regression test result with explanatory messages
- The test checks that optimcon() accepts a minimal one-spin Hilbert-space
- control problem, and that grape_hilb() fidelity and gradient agree with
- independent matrix exponentiation and finite differences.
- Announce the test target
- State the optimal-control target of the test
- Ensure that a parallel pool is available for the ensemble loop
- Build a minimal Hilbert-space Spinach object for setup and GRAPE
- Define one-spin operators and a two-step timing grid
- Configure a minimal optimal-control problem

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `optimcon()`, `grape_hilb()`, `gcp()`, `parpool()`, `local_spin_system()`, `pauli()`, `test_true()`, `test_close()`, `rmfield()`, `strcmp()`, `contains()`, `waveform()`, `pulse_dt()`, `hdot()`, `wf_plus()`.
