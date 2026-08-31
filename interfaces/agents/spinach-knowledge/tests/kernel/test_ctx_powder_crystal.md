# tests/kernel/test_ctx_powder_crystal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_powder_crystal.m`
- Signature: `result=test_ctx_powder_crystal()`
- Total lines: 58

## Purpose

Tests static powder and crystal contexts at one orientation. Syntax: result=test_ctx_powder_crystal()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Powder and crystal single orientation\n')`.
- Lines 20-23: State the static-context target of the test; implemented by `result=new_test_result('kernel/ctx_powder_crystal', 'Powder and crystal single orientation', 'powder() with the single_crystal grid must match crystal().')`.
- Lines 25-26: Build a one-spin anisotropic Liouville-space system; implemented by `sys.magnet=14.1`.
- Lines 35-36: Set up a short static acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 48-49: Run both production static contexts; implemented by `fid_powder=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 52-54: Check that one-point powder averaging is a crystal calculation; implemented by `result=test_close(result,'single orientation powder',fid_powder,fid_crystal,1e-12,1e-12, 'the single_crystal powder grid has unit weight at zero Euler angles')`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/ctx_powder_crystal', 'Powder and crystal single orientation', 'powder() with the single_crystal grid must match crystal().')`.
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
- Lines 42: computes `parameters.npoints` using `parameters.npoints=4`.
- Lines 43: computes `parameters.grid` using `parameters.grid='single_crystal'`.
- Lines 44: computes `parameters.orientation` using `parameters.orientation=[0 0 0]`.

## Outputs

- result -regression test result with explanatory messages
- The test uses the single_crystal grid so that powder() and crystal()
- represent the same Euler orientation of the same anisotropic one-spin
- Hamiltonian.

## Implementation structure

- Tests static powder and crystal contexts at one orientation. Syntax:
- result=test_ctx_powder_crystal()
- result -regression test result with explanatory messages
- The test uses the single_crystal grid so that powder() and crystal()
- represent the same Euler orientation of the same anisotropic one-spin
- Hamiltonian.
- Announce the test target
- State the static-context target of the test
- Build a one-spin anisotropic Liouville-space system
- Set up a short static acquisition
- Run both production static contexts
- Check that one-point powder averaging is a crystal calculation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `powder()`, `crystal()`, `test_spin_system()`, `state()`, `test_close()`.
