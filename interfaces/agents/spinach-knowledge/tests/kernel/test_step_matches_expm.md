# tests/kernel/test_step_matches_expm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_step_matches_expm.m`
- Signature: `result=test_step_matches_expm()`
- Total lines: 48

## Purpose

Tests Hilbert-space propagation against matrix exponentiation. Syntax: result=test_step_matches_expm()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Hilbert propagation against expm\n')`.
- Lines 19-22: State the propagation target of the test; implemented by `result=new_test_result('kernel/step_matches_expm', 'Hilbert propagation against expm', 'step() must reproduce unitary density-matrix propagation.')`.
- Lines 24-25: Build a one-proton Hilbert-space spin system; implemented by `sys.magnet=0`.
- Lines 32-33: Define a Hamiltonian and an initial density matrix; implemented by `S=pauli(2)`.
- Lines 38-39: Build the independent exact propagator; implemented by `P=expm(-1i*H*dt)`.
- Lines 43-45: Check exact finite-dimensional propagation; implemented by `result=test_close(result,'step versus expm',rho_obs,rho_ref,1e-13,1e-13, 'finite Hilbert-space propagation is exactly unitary')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/step_matches_expm', 'Hilbert propagation against expm', 'step() must reproduce unitary density-matrix propagation.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 33: computes `S` using `S=pauli(2)`.
- Lines 34: computes `H` using `H=2*pi*123*S.z`.
- Lines 35: computes `rho` using `rho=S.x+0.25*S.y`.
- Lines 36: computes `dt` using `dt=2.5e-3`.
- Lines 39: computes `P` using `P=expm(-1i*H*dt)`.
- Lines 40: computes `rho_ref` using `rho_ref=P*rho*P'`.
- Lines 41: computes `rho_obs` using `rho_obs=step(spin_system,H,rho,dt)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks the Spinach sign convention for density-matrix evolution:
- rho(t)=exp(-iHt) rho(0) exp(+iHt).

## Implementation structure

- Tests Hilbert-space propagation against matrix exponentiation. Syntax:
- result=test_step_matches_expm()
- result -regression test result with explanatory messages
- The test checks the Spinach sign convention for density-matrix evolution:
- rho(t)=exp(-iHt) rho(0) exp(+iHt).
- Announce the test target
- State the propagation target of the test
- Build a one-proton Hilbert-space spin system
- Define a Hamiltonian and an initial density matrix
- Build the independent exact propagator
- Check exact finite-dimensional propagation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `step()`, `test_spin_system()`, `pauli()`, `test_close()`.
