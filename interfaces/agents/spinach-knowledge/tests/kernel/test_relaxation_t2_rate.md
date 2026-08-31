# tests/kernel/test_relaxation_t2_rate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_relaxation_t2_rate.m`
- Signature: `result=test_relaxation_t2_rate()`
- Total lines: 48

## Purpose

Tests phenomenological T2 relaxation rate. Syntax: result=test_relaxation_t2_rate()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Phenomenological T2 decay rate\n')`.
- Lines 19-22: State the relaxation target of the test; implemented by `result=new_test_result('kernel/relaxation_t2_rate', 'Phenomenological T2 decay rate', 'the t1_t2 model must assign transverse magnetisation the generator eigenvalue -R2.…`.
- Lines 24-25: Build a one-spin relaxation system; implemented by `sys.magnet=14.1`.
- Lines 38-39: Build relaxation superoperator and transverse state; implemented by `R=relaxation(spin_system)`.
- Lines 43-45: Check that L+ is an eigenstate with the negative R2 generator eigenvalue; implemented by `result=test_close(result,'R2 eigenvalue on L+',Rrho,-7.0*rho,1e-12,1e-12, 'decay generators have negative eigenvalues, with magnitude R2')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/relaxation_t2_rate', 'Phenomenological T2 decay rate', 'the t1_t2 model must assign transverse magnetisation the generator eigenvalue -R2.…`.
- Lines 25: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 29: computes `inter.r1_rates` using `inter.r1_rates={2.0}`.
- Lines 30: computes `inter.r2_rates` using `inter.r2_rates={7.0}`.
- Lines 31: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 32: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 33: computes `inter.temperature` using `inter.temperature=298`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 39: computes `R` using `R=relaxation(spin_system)`.
- Lines 40: computes `rho` using `rho=state(spin_system,'L+','1H')`.
- Lines 41: computes `Rrho` using `Rrho=R*rho`.

## Outputs

- result -regression test result with explanatory messages
- The test checks that the t1_t2 relaxation model gives transverse L+ order
- the negative generator eigenvalue corresponding to the specified R2 rate.

## Implementation structure

- Tests phenomenological T2 relaxation rate. Syntax:
- result=test_relaxation_t2_rate()
- result -regression test result with explanatory messages
- The test checks that the t1_t2 relaxation model gives transverse L+ order
- the negative generator eigenvalue corresponding to the specified R2 rate.
- Announce the test target
- State the relaxation target of the test
- Build a one-spin relaxation system
- Build relaxation superoperator and transverse state
- Check that L+ is an eigenstate with the negative R2 generator eigenvalue

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `relaxation()`, `state()`, `test_close()`.
