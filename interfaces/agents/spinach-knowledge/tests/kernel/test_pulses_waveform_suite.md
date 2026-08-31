# tests/kernel/test_pulses_waveform_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_pulses_waveform_suite.m`
- Signature: `result=test_pulses_waveform_suite()`
- Total lines: 126

## Purpose

Tests deterministic pulse waveform generators. Syntax: result=test_pulses_waveform_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Pulse waveform generators\n')`.
- Lines 19-22: State the waveform target of the test; implemented by `result=new_test_result('kernel/pulses_waveform_suite', 'Pulse waveform generators', 'pulse waveform utilities must reproduce their analytic formulae and tabulated period…`.
- Lines 24-25: Check sawtooth and triangular wave formulae at simple fractions of the period; implemented by `amp=2`.
- Lines 35-36: Check Uhrig delay formula for three pulses over a unit interval; implemented by `T=1`.
- Lines 45-47: Check periodic phase tables and wrap-around indexing; implemented by `result=test_close(result,'pmlg5 first phase',pmlg5(1),(pi/180)*339.22,1e-14,1e-14, 'PMLG5 index one is the first tabulated phase')`.
- Lines 55-57: Check simple analytic pulse envelopes; implemented by `result=test_close(result,'pulse_shape rectangular',pulse_shape('rectangular',4),ones(1,4),1e-15,1e-15, 'the rectangular envelope is constant over all pulse slices')`.
- Lines 61-62: Check reading of a distributed rectangular Bruker pulse file; implemented by `[A,phi,Cx,Cy,scaling_factor]=read_wave('rectangular_1000.pk',4)`.
- Lines 74-75: Check Veshtort-Griffin duration scaling without duplicating the coefficient table; implemented by `vg_short=vg_pulse('E0A',7,1)`.
- Lines 80-81: Check WURST chirp formula on a grid whose samples are simple fractions; implemented by `[Cx,Cy,durs,ints,amps,phis,frqs]=chirp_pulse(5,1,4,2,'wurst')`.
- Lines 102-103: Check hyperbolic secant pulse formulae; implemented by `peak_amp=3`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/pulses_waveform_suite', 'Pulse waveform generators', 'pulse waveform utilities must reproduce their analytic formulae and tabulated period…`.
- Lines 25: computes `amp` using `amp=2`.
- Lines 26: computes `freq` using `freq=1`.
- Lines 27: computes `time_grid` using `time_grid=[0 0.25 0.5 0.75 1.0]`.
- Lines 28: computes `saw_ref` using `saw_ref=[-2 -1 0 1 -2]`.
- Lines 29: computes `tri_ref` using `tri_ref=abs(saw_ref)`.
- Lines 36: computes `T` using `T=1`.
- Lines 37: computes `N` using `N=3`.
- Lines 38: computes `pos` using `pos=T*(sin(pi*(1:N)/(2*N+2)).^2-0.5)`.
- Lines 39: computes `delays` using `delays=diff(pos)`.
- Lines 40: computes `chunk` using `chunk=(T-sum(delays))/2`.
- Lines 41: computes `udd_ref` using `udd_ref=[chunk delays chunk]`.
- Lines 62: computes `[A,phi,Cx,Cy,scaling_factor]` using `[A,phi,Cx,Cy,scaling_factor]=read_wave('rectangular_1000.pk',4)`.
- Lines 75: computes `vg_short` using `vg_short=vg_pulse('E0A',7,1)`.
- Lines 76: computes `vg_long` using `vg_long=vg_pulse('E0A',7,2)`.
- Lines 81: computes `[Cx,Cy,durs,ints,amps,phis,frqs]` using `[Cx,Cy,durs,ints,amps,phis,frqs]=chirp_pulse(5,1,4,2,'wurst')`.
- Lines 83: computes `amps_ref` using `amps_ref=2*pi*sqrt(4)*(1-abs(sin(pi*time_grid).^2))`.
- Lines 84: computes `phis_ref` using `phis_ref=pi*4*(time_grid.^2)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks analytic waveform formulas, periodic phase tables,
- JCAMP pulse-file reading, chirp construction, and sech pulse coordinates.

## Implementation structure

- Tests deterministic pulse waveform generators. Syntax:
- result=test_pulses_waveform_suite()
- result -regression test result with explanatory messages
- The test checks analytic waveform formulas, periodic phase tables,
- JCAMP pulse-file reading, chirp construction, and sech pulse coordinates.
- Announce the test target
- State the waveform target of the test
- Check sawtooth and triangular wave formulae at simple fractions of the period
- Check Uhrig delay formula for three pulses over a unit interval
- Check periodic phase tables and wrap-around indexing
- Check simple analytic pulse envelopes
- Check reading of a distributed rectangular Bruker pulse file

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `sawtooth()`, `triwave()`, `diff()`, `uhrig_times()`, `pmlg5()`, `spinal()`, `pulse_shape()`, `read_wave()`, `vg_pulse()`, `chirp_pulse()`, `polar2cartesian()`, `sech_pulse()`, `sech()`, `cosh()`.
