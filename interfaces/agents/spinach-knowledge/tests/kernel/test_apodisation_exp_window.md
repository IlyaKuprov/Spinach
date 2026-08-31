# tests/kernel/test_apodisation_exp_window.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_apodisation_exp_window.m`
- Signature: `result=test_apodisation_exp_window()`
- Total lines: 38

## Purpose

Tests exponential FID apodisation. Syntax: result=test_apodisation_exp_window()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Exponential FID apodisation\n')`.
- Lines 19-22: State the processing target of the test; implemented by `result=new_test_result('kernel/apodisation_exp_window', 'Exponential FID apodisation', 'exponential apodisation must multiply by exp(-k*x) and halve the first point.')`.
- Lines 24-25: Build a minimal reporting object and a constant FID; implemented by `spin_system.sys.output='hush'`.
- Lines 28-29: Apply an exponential window; implemented by `fid_obs=apodisation(spin_system,fid,{{'exp',1}})`.
- Lines 33-35: Check the explicit window; implemented by `result=test_close(result,'exp window and first-point half',fid_obs,fid_ref,1e-15,1e-15, 'the first point is halved, then multiplied by exp(-x)')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/apodisation_exp_window', 'Exponential FID apodisation', 'exponential apodisation must multiply by exp(-k*x) and halve the first point.')`.
- Lines 25: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 26: computes `fid` using `fid=ones(4,1)`.
- Lines 29: computes `fid_obs` using `fid_obs=apodisation(spin_system,fid,{{'exp',1}})`.
- Lines 30: computes `fid_ref` using `fid_ref=exp(-linspace(0,1,4)).'`.
- Lines 31: computes `fid_ref(1)` using `fid_ref(1)=fid_ref(1)/2`.

## Outputs

- result -regression test result with explanatory messages
- The test checks the explicit exponential window and the NMR Fourier
- convention that halves the first point of each active FID dimension.

## Implementation structure

- Tests exponential FID apodisation. Syntax:
- result=test_apodisation_exp_window()
- result -regression test result with explanatory messages
- The test checks the explicit exponential window and the NMR Fourier
- convention that halves the first point of each active FID dimension.
- Announce the test target
- State the processing target of the test
- Build a minimal reporting object and a constant FID
- Apply an exponential window
- Check the explicit window

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `apodisation()`, `fid_ref()`, `test_close()`.
