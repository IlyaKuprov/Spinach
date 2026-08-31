# tests/kernel/test_liquid_single_spin_fid.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_liquid_single_spin_fid.m`
- Signature: `result=test_liquid_single_spin_fid()`
- Total lines: 52

## Purpose

Tests a one-spin liquid-state free induction decay. Syntax: result=test_liquid_single_spin_fid()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Single-spin liquid-state FID\n')`.
- Lines 19-22: State the NMR target of the test; implemented by `result=new_test_result('kernel/liquid_single_spin_fid', 'Single-spin liquid-state FID', 'a zero-offset isolated spin has a constant free induction decay.')`.
- Lines 24-25: Build a one-spin Liouville-space system; implemented by `sys.magnet=14.1`.
- Lines 32-33: Set up a zero-offset acquisition; implemented by `parameters.spins={'1H'}`.
- Lines 44-45: Simulate the FID; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 47-49: The signal is constant and equal to its first point; implemented by `result=test_close(result,'constant zero-offset FID',fid,fid(1)*ones(size(fid)),1e-12,1e-12, 'without precession or relaxation the detected coherence is time-independent')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/liquid_single_spin_fid', 'Single-spin liquid-state FID', 'a zero-offset isolated spin has a constant free induction decay.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 33: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 34: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 35: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 36: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 37: computes `parameters.offset` using `parameters.offset=0`.
- Lines 38: computes `parameters.sweep` using `parameters.sweep=1000`.
- Lines 39: computes `parameters.npoints` using `parameters.npoints=8`.
- Lines 40: computes `parameters.zerofill` using `parameters.zerofill=8`.
- Lines 41: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 42: computes `parameters.invert_axis` using `parameters.invert_axis=0`.
- Lines 45: computes `fid` using `fid=liquid(spin_system,@acquire,parameters,'nmr')`.

## Outputs

- result -regression test result with explanatory messages
- The test simulates a zero-offset one-spin FID. With no Hamiltonian and no
- relaxation, transverse magnetisation is constant in time.

## Implementation structure

- Tests a one-spin liquid-state free induction decay. Syntax:
- result=test_liquid_single_spin_fid()
- result -regression test result with explanatory messages
- The test simulates a zero-offset one-spin FID. With no Hamiltonian and no
- relaxation, transverse magnetisation is constant in time.
- Announce the test target
- State the NMR target of the test
- Build a one-spin Liouville-space system
- Set up a zero-offset acquisition
- Simulate the FID
- The signal is constant and equal to its first point

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `state()`, `liquid()`, `test_close()`, `fid()`.
