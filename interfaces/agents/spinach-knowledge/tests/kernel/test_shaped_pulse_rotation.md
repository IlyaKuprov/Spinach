# tests/kernel/test_shaped_pulse_rotation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_shaped_pulse_rotation.m`
- Signature: `result=test_shaped_pulse_rotation()`
- Total lines: 49

## Purpose

Tests a one-slice Cartesian shaped pulse. Syntax: result=test_shaped_pulse_rotation()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Piecewise Cartesian pulse rotation\n')`.
- Lines 19-22: State the pulse target of the test; implemented by `result=new_test_result('kernel/shaped_pulse_rotation', 'Piecewise Cartesian pulse rotation', 'a rectangular Cartesian shaped pulse must reproduce the hard-pulse limit.')`.
- Lines 24-25: Build a one-proton Hilbert-space spin system; implemented by `sys.magnet=0`.
- Lines 32-33: Define drift, controls, and one pi pulse slice; implemented by `Lx=operator(spin_system,'Lx',1)`.
- Lines 41-42: Apply the shaped pulse; implemented by `rho_obs=shaped_pulse_xy(spin_system,drift,controls,amplitudes,slice_durs,Lz,'expm-pwc')`.
- Lines 44-46: Check the physical rotation result; implemented by `result=test_close(result,'rectangular X pi pulse',rho_obs,-Lz,1e-14,1e-14, 'amplitude times duration is the flip angle in radians')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/shaped_pulse_rotation', 'Piecewise Cartesian pulse rotation', 'a rectangular Cartesian shaped pulse must reproduce the hard-pulse limit.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 33: computes `Lx` using `Lx=operator(spin_system,'Lx',1)`.
- Lines 34: computes `Ly` using `Ly=operator(spin_system,'Ly',1)`.
- Lines 35: computes `Lz` using `Lz=state(spin_system,'Lz',1)`.
- Lines 36: computes `drift` using `drift=0*Lx`.
- Lines 37: computes `controls` using `controls={Lx,Ly}`.
- Lines 38: computes `amplitudes` using `amplitudes={1,0}`.
- Lines 39: computes `slice_durs` using `slice_durs=pi`.
- Lines 42: computes `rho_obs` using `rho_obs=shaped_pulse_xy(spin_system,drift,controls,amplitudes,slice_durs,Lz,'expm-pwc')`.

## Outputs

- result -regression test result with explanatory messages
- The test applies one rectangular X pulse slice with amplitude 1 rad/s and
- duration pi seconds; the net flip angle is pi, so Lz must invert.

## Implementation structure

- Tests a one-slice Cartesian shaped pulse. Syntax:
- result=test_shaped_pulse_rotation()
- result -regression test result with explanatory messages
- The test applies one rectangular X pulse slice with amplitude 1 rad/s and
- duration pi seconds; the net flip angle is pi, so Lz must invert.
- Announce the test target
- State the pulse target of the test
- Build a one-proton Hilbert-space spin system
- Define drift, controls, and one pi pulse slice
- Apply the shaped pulse
- Check the physical rotation result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `operator()`, `state()`, `shaped_pulse_xy()`, `test_close()`.
