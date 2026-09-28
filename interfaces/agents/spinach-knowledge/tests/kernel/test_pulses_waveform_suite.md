# tests/kernel/test_pulses_waveform_suite.m

- Signature: `result=test_pulses_waveform_suite()`

Regression-tests waveform utilities against analytic values. `result` contains test messages.

- For `amp=2`, `freq=1`, and `t=[0,.25,.5,.75,1]`, checks sawtooth `amp*(2*freq*mod(t,1/freq)-1)` and triangle-wave absolute value.
- For `T=1`, `N=3`, checks Uhrig positions `pos=T*(sin(pi*(1:N)/(2*N+2)).^2-0.5)`, delays `diff(pos)`, and equal end chunks `(T-sum(delays))/2`.
- Checks PMLG5/SPINAL first phases `339.22°`/`10°` and periods `20`/`64`.
- Checks rectangular and sinc3 envelopes; `rectangular_1000.pk` amplitude, phase, Cartesian controls and scaling; and VG `E0A` inverse-duration scaling.
- Checks WURST chirp (5 points, duration 1, bandwidth 4, exponent 2) and sech pulse (3, 2, 5, 2, 5): grids, amplitudes, phases and Cartesian controls.
