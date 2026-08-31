# tests/kernel/test_slr_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_slr_pulse.m`
- Signature: `result=test_slr_pulse()`
- Total lines: 191

## Purpose

Tests Shinnar-Le Roux selective excitation pulse design. Syntax: result=test_slr_pulse()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `throws_with()`, `ck_profile()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Shinnar-Le Roux pulse design\n')`.
- Lines 20-23: State the selective pulse target of the test; implemented by `result=new_test_result('kernel/slr_pulse', 'Shinnar-Le Roux pulse design', 'an SLR waveform must produce the requested flip and a selective excitation profile.')`.
- Lines 25-26: Define a representative selective excitation design; implemented by `npts=64`.
- Lines 33-34: Generate the production waveform; implemented by `[Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)`.
- Lines 36-41: Check the output dimensions and finiteness; implemented by `result=test_true(result,'SLR output dimensions', isequal(size(Cx),[1 npts])&&isequal(size(Cy),[1 npts])&& isequal(size(durs),[1 npts])&&isequal(size(amps),[1 npts])&& is…`.
- Lines 48-50: Check Cartesian and polar coordinate consistency; implemented by `result=test_close(result,'SLR polar X consistency',amps.*cos(phis),Cx,1e-10,1e-14, 'polar amplitude and phase must reproduce the Cartesian X control')`.
- Lines 54-55: Build independent spin-half matrices; implemented by `Lx=[0 1;1 0]/2`.
- Lines 59-60: Propagate the waveform directly at zero offset; implemented by `rho=Lz`.
- Lines 66-68: Check the requested flip and Spinach phase convention; implemented by `result=test_close(result,'SLR direct centre flip',rho,-Ly,1e-11,1e-11, 'a positive-X pi/2 control under exp(-1i*H*t) rotates Lz to -Ly')`.
- Lines 70-71: Check an interior value of the documented flip-angle interval; implemented by `check_flip=pi/6`.
- Lines 83-84: Build a dense independent Cayley-Klein frequency sweep; implemented by `freq_grid=linspace(-0.5,0.5,16385)`.
- Lines 87-88: Reconstruct the documented transition edges independently; implemented by `beta_pass=sqrt(pass_rip/2)`.
- Lines 99-101: Check unitarity and the selective excitation profile; implemented by `result=test_true(result,'SLR Cayley-Klein unitarity',unit_err<2e-12, 'every frequency response must remain unitary throughout propagation')`.
- Lines 109-110: Check that a smaller flip is designed before the nonlinear inverse transform; implemented by `[small_trans,small_long,small_unit]=ck_profile(small_x,small_y,small_durs,freq_grid)`.
- Lines 122-123: Build a one-proton Hilbert-space spin system; implemented by `sys.magnet=0`.
- Lines 130-131: Apply the generated controls through the production shaped-pulse path; implemented by `Lx_op=operator(spin_system,'Lx',1)`.
- Lines 138-140: Check production-path sign and units; implemented by `result=test_close(result,'SLR shaped-pulse centre flip',rho_obs,-Ly_state,1e-10,1e-10, 'the generated rad/s controls and second durations must produce the requested flip…`.
- Lines 142-145: Check representative consistency failures; implemented by `result=test_true(result,'SLR odd sample count rejected', throws_with(@()slr_pulse(63,dur,tbw,flip_angle,pass_rip,stop_rip),'even'), 'the linear-phase design requires an…`.

### Control flow inferred from the code

- Line 61: `for` loop over `slice=1:npts`.
- Line 74: `for` loop over `slice=1:npts`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/slr_pulse', 'Shinnar-Le Roux pulse design', 'an SLR waveform must produce the requested flip and a selective excitation profile.')`.
- Lines 26: computes `npts` using `npts=64`.
- Lines 27: computes `dur` using `dur=4e-3`.
- Lines 28: computes `tbw` using `tbw=4`.
- Lines 29: computes `flip_angle` using `flip_angle=pi/2`.
- Lines 30: computes `pass_rip` using `pass_rip=0.01`.
- Lines 31: computes `stop_rip` using `stop_rip=0.01`.
- Lines 34: computes `[Cx,Cy,durs,amps,phis]` using `[Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)`.
- Lines 55: computes `Lx` using `Lx=[0 1;1 0]/2`.
- Lines 56: computes `Ly` using `Ly=[0 -1i;1i 0]/2`.
- Lines 57: computes `Lz` using `Lz=[1 0;0 -1]/2`.
- Lines 60: computes `rho` using `rho=Lz`.
- Lines 62: computes `U` using `U=expm(-1i*(Cx(slice)*Lx+Cy(slice)*Ly)*durs(slice))`.
- Lines 71: computes `check_flip` using `check_flip=pi/6`.
- Lines 72: computes `[small_x,small_y,small_durs]` using `[small_x,small_y,small_durs]=slr_pulse(npts,dur,tbw,check_flip,pass_rip,stop_rip)`.
- Lines 78: computes `flip_ref` using `flip_ref=cos(check_flip)*Lz-sin(check_flip)*Ly`.
- Lines 84: computes `freq_grid` using `freq_grid=linspace(-0.5,0.5,16385)`.
- Lines 85: computes `[transverse,longitudinal,unit_err]` using `[transverse,longitudinal,unit_err]=ck_profile(Cx,Cy,durs,freq_grid)`.

### Local helper functions

- Line 162: `throws_with()` — `function verdict=throws_with(fun_handle,message_text)`. Evaluate a waveform with an independent Cayley-Klein propagator
  - Representative operation: `try`.
  - Representative operation: `fun_handle()`.
- Line 172: `ck_profile()` — `function [transverse,longitudinal,unit_err]=ck_profile(Cx,Cy,durs,freq_grid)`.
  - Representative operation: `offsets=2*pi*freq_grid/durs(1)`.
  - Representative operation: `alpha=ones(size(freq_grid))`.

## Outputs

- result -regression test result with explanatory messages
- The test checks waveform units and shape, independent two-level
- propagation, excitation profile selectivity, production-path shaped
- pulse propagation, and representative input validation failures.

## Implementation structure

- Tests Shinnar-Le Roux selective excitation pulse design. Syntax:
- result=test_slr_pulse()
- result -regression test result with explanatory messages
- The test checks waveform units and shape, independent two-level
- propagation, excitation profile selectivity, production-path shaped
- pulse propagation, and representative input validation failures.
- Announce the test target
- State the selective pulse target of the test
- Define a representative selective excitation design
- Generate the production waveform
- Check the output dimensions and finiteness
- Check Cartesian and polar coordinate consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `slr_pulse()`, `test_true()`, `isequal()`, `all()`, `test_close()`, `durs()`, `small_x()`, `small_y()`, `small_durs()`, `num2str()`, `ck_profile()`, `log10()`, `transverse()`, `longitudinal()`, `small_trans()`.
