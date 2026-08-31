# tests/kernel/test_dynamic_relaxation_models_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_relaxation_models_suite.m`
- Signature: `result=test_dynamic_relaxation_models_suite()`
- Total lines: 160

## Purpose

Tests dynamic relaxation model helper paths. Syntax: result=test_dynamic_relaxation_models_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Dynamic relaxation model helpers\n')`.
- Lines 19-22: State the relaxation-model target of the test; implemented by `result=new_test_result('kernel/dynamic_relaxation_models_suite', 'Dynamic relaxation model helpers', 'relaxation helpers must assign mathematically expected rates on tin…`.
- Lines 24-25: Check anisotropic tensor rates in the extended T1/T2 model; implemented by `sys.magnet=14.1`.
- Lines 45-46: Check function-handle rates in the extended T1/T2 model; implemented by `inter.r1_rates={@(alp,bet,gam)2+0*alp+0*bet+0*gam}`.
- Lines 57-58: Check non-selective damping in Liouville space; implemented by `inter_d.zeeman.scalar={0}`.
- Lines 73-74: Check one-spin Lindblad relaxation rates; implemented by `inter_l.zeeman.scalar={0}`.
- Lines 93-94: Check scalar Redfield integral against the closed zero-H0 reference; implemented by `spin_l.tols.rlx_integration=1e-5`.
- Lines 103-104: Build a Redfield-ready spin system for correlation-function checks; implemented by `sys_r.magnet=14.1`.
- Lines 116-117: Check isotropic rank-two correlation-function weights and rates; implemented by `[weights,rates,states]=corrfun(spin_r,2,3,2,3,2)`.
- Lines 125-126: Check axial rotational diffusion rate formula; implemented by `spin_r.rlx.tau_c={[2e-9 5e-9]}`.
- Lines 136-137: Check rhombic rotational diffusion rate formula; implemented by `spin_r.rlx.tau_c={[1e-9 2e-9 3e-9]}`.
- Lines 148-149: Exercise the serial Redfield include on a tiny anisotropic one-spin system; implemented by `inter_r.zeeman.matrix={diag([1 2 -3])}`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_relaxation_models_suite', 'Dynamic relaxation model helpers', 'relaxation helpers must assign mathematically expected rates on tin…`.
- Lines 25: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 29: computes `inter.r1_rates` using `inter.r1_rates={diag([1 3 5])}`.
- Lines 30: computes `inter.r2_rates` using `inter.r2_rates={diag([2 4 6])}`.
- Lines 31: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 32: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 33: computes `inter.temperature` using `inter.temperature=300`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_t` using `spin_t=test_spin_system(sys,inter,bas)`.
- Lines 37: computes `[r1_op,r2_op]` using `[r1_op,r2_op]=rlx_t1_t2(spin_t,[0 0 0])`.
- Lines 38: computes `rho_z` using `rho_z=state(spin_t,'Lz','1H')`.
- Lines 39: computes `rho_p` using `rho_p=state(spin_t,'L+','1H')`.
- Lines 48: computes `spin_f` using `spin_f=test_spin_system(sys,inter,bas)`.
- Lines 58: computes `inter_d.zeeman.scalar` using `inter_d.zeeman.scalar={0}`.

## Outputs

- result -regression test result with explanatory messages
- The test checks anisotropic and functional T1/T2 rates, damping, Lindblad,
- scalar Redfield, correlation functions, and the serial Redfield include.

## Implementation structure

- Tests dynamic relaxation model helper paths. Syntax:
- result=test_dynamic_relaxation_models_suite()
- result -regression test result with explanatory messages
- The test checks anisotropic and functional T1/T2 rates, damping, Lindblad,
- scalar Redfield, correlation functions, and the serial Redfield include.
- Announce the test target
- State the relaxation-model target of the test
- Check anisotropic tensor rates in the extended T1/T2 model
- Check function-handle rates in the extended T1/T2 model
- Check non-selective damping in Liouville space
- Check one-spin Lindblad relaxation rates
- Check scalar Redfield integral against the closed zero-H0 reference

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `rlx_t1_t2()`, `state()`, `test_close()`, `relaxation()`, `unit_state()`, `rlx_scalar()`, `speye()`, `corrfun()`, `nnz()`, `test_true()`.
