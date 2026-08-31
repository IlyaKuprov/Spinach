# tests/kernel/test_ctx_singlerot_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_singlerot_acquire.m`
- Signature: `result=test_ctx_singlerot_acquire()`
- Total lines: 68

## Purpose

Tests the single-rotor context with acquire(). Syntax: result=test_ctx_singlerot_acquire()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Single-rotor acquire path\n')`.
- Lines 20-23: State the single-rotor target of the test; implemented by `result=new_test_result('kernel/ctx_singlerot_acquire', 'Single-rotor acquire path', 'singlerot() must project states into rotor space and run acquire().')`.
- Lines 25-26: Build a one-spin anisotropic Liouville-space system; implemented by `sys.magnet=14.1`.
- Lines 35-36: Set up a tiny MAS acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 50-51: Run the production single-rotor context; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 53-55: Check the number of acquired points; implemented by `result=test_close(result,'singlerot FID length',numel(fid),parameters.npoints,0,0, 'acquire() should return one point for each requested time sample')`.
- Lines 57-58: Check that the zero-time signal survives state projection; implemented by `fid_zero=parameters.coil'*parameters.rho0`.
- Lines 62-64: Check that the acquired trace is finite; implemented by `result=test_true(result,'singlerot finite FID',all(isfinite(fid(:))), 'short MAS propagation should not produce NaN or Inf values')`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/ctx_singlerot_acquire', 'Single-rotor acquire path', 'singlerot() must project states into rotor space and run acquire().')`.
- Lines 26: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 27: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 28: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[-2 -2 4]}`.
- Lines 29: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0]}`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `bas.projections` using `bas.projections=+1`.
- Lines 33: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 37: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 38: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 39: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 40: computes `parameters.offset` using `parameters.offset=0`.
- Lines 41: computes `parameters.sweep` using `parameters.sweep=2000`.
- Lines 42: computes `parameters.npoints` using `parameters.npoints=3`.
- Lines 43: computes `parameters.rate` using `parameters.rate=1000`.
- Lines 44: computes `parameters.axis` using `parameters.axis=[1 1 1]`.

## Outputs

- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin MAS calculation through
- singlerot() and checks the returned time-domain trace for basic
- physical and dimensional invariants.

## Implementation structure

- Tests the single-rotor context with acquire(). Syntax:
- result=test_ctx_singlerot_acquire()
- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin MAS calculation through
- singlerot() and checks the returned time-domain trace for basic
- physical and dimensional invariants.
- Announce the test target
- State the single-rotor target of the test
- Build a one-spin anisotropic Liouville-space system
- Set up a tiny MAS acquisition
- Run the production single-rotor context
- Check the number of acquired points

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `singlerot()`, `acquire()`, `test_spin_system()`, `state()`, `test_close()`, `fid()`, `test_true()`, `all()`.
