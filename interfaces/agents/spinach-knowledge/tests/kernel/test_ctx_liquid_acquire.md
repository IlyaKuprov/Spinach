# tests/kernel/test_ctx_liquid_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_liquid_acquire.m`
- Signature: `result=test_ctx_liquid_acquire()`
- Total lines: 62

## Purpose

Tests the liquid context against the direct acquire path. Syntax: result=test_ctx_liquid_acquire()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Liquid context acquire path\n')`.
- Lines 20-23: State the liquid-context target of the test; implemented by `result=new_test_result('kernel/ctx_liquid_acquire', 'Liquid context acquire path', 'liquid() must pass the offset Liouvillian to acquire().')`.
- Lines 25-26: Build a one-spin Liouville-space system; implemented by `sys.magnet=14.1`.
- Lines 34-35: Set up a short offset acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 45-46: Run the production liquid context; implemented by `fid_ctx=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 48-49: Build the same direct acquire input path; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 56-58: Check that the context passes through the direct dynamics; implemented by `result=test_close(result,'liquid context FID',fid_ctx,fid_ref,1e-12,1e-12, 'liquid() should reproduce direct acquire() for the same offset Liouvillian')`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/ctx_liquid_acquire', 'Liquid context acquire path', 'liquid() must pass the offset Liouvillian to acquire().')`.
- Lines 26: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 27: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 28: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `bas.projections` using `bas.projections=+1`.
- Lines 32: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 35: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 36: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 37: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 38: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 39: computes `parameters.offset` using `parameters.offset=25`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=1000`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=4`.
- Lines 42: computes `parameters.needs` using `parameters.needs={}`.
- Lines 43: computes `parameters.rframes` using `parameters.rframes={}`.
- Lines 46: computes `fid_ctx` using `fid_ctx=liquid(spin_system,@acquire,parameters,'nmr')`.

## Outputs

- result -regression test result with explanatory messages
- The test runs a one-spin offset FID through liquid() and compares it
- with the same Hamiltonian, relaxation, and kinetics objects passed
- directly to acquire().

## Implementation structure

- Tests the liquid context against the direct acquire path. Syntax:
- result=test_ctx_liquid_acquire()
- result -regression test result with explanatory messages
- The test runs a one-spin offset FID through liquid() and compares it
- with the same Hamiltonian, relaxation, and kinetics objects passed
- directly to acquire().
- Announce the test target
- State the liquid-context target of the test
- Build a one-spin Liouville-space system
- Set up a short offset acquisition
- Run the production liquid context
- Build the same direct acquire input path

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `liquid()`, `acquire()`, `test_spin_system()`, `state()`, `assume()`, `hamiltonian()`, `relaxation()`, `kinetics()`, `frqoffset()`, `test_close()`.
