# tests/kernel/test_step_zero_time.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_step_zero_time.m`
- Signature: `result=test_step_zero_time()`
- Total lines: 43

## Purpose

Tests zero-duration propagation. Syntax: result=test_step_zero_time()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Zero-duration propagation identity\n')`.
- Lines 19-22: State the propagation target of the test; implemented by `result=new_test_result('kernel/step_zero_time', 'Zero-duration propagation identity', 'a propagator over zero time is the identity map.')`.
- Lines 24-25: Build a one-proton Hilbert-space spin system; implemented by `sys.magnet=0`.
- Lines 32-33: Propagate an arbitrary Hermitian density matrix for zero time; implemented by `S=pauli(2)`.
- Lines 38-40: Check the identity limit; implemented by `result=test_close(result,'rho(t=0)=rho(0)',rho_obs,rho,1e-15,1e-15, 'zero-duration evolution cannot change the state')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/step_zero_time', 'Zero-duration propagation identity', 'a propagator over zero time is the identity map.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 33: computes `S` using `S=pauli(2)`.
- Lines 34: computes `rho` using `rho=S.x+2*S.z`.
- Lines 35: computes `H` using `H=3*S.x+5*S.z`.
- Lines 36: computes `rho_obs` using `rho_obs=step(spin_system,H,rho,0)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks the identity limit of the propagator: a zero time step
- must leave the density matrix exactly unchanged.

## Implementation structure

- Tests zero-duration propagation. Syntax:
- result=test_step_zero_time()
- result -regression test result with explanatory messages
- The test checks the identity limit of the propagator: a zero time step
- must leave the density matrix exactly unchanged.
- Announce the test target
- State the propagation target of the test
- Build a one-proton Hilbert-space spin system
- Propagate an arbitrary Hermitian density matrix for zero time
- Check the identity limit

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `pauli()`, `step()`, `test_close()`, `rho()`.
