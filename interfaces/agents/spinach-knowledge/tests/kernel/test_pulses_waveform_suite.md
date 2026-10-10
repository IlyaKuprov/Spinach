# tests/kernel/test_pulses_waveform_suite.m

## Purpose

Regression test for the deterministic pulse waveform generators in Spinach. The suite verifies that pulse waveform utilities reproduce their analytic formulae and tabulated periodic sequences, covering sawtooth and triangular wave formulae, the Uhrig delay formula, periodic phase tables, analytic pulse envelopes, JCAMP pulse-file reading, Veshtort-Griffin duration scaling, WURST chirp construction, and hyperbolic secant pulse coordinates.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_pulses_waveform_suite.m>

## Behaviour

The function announces the test target with `fprintf('TESTING: Pulse waveform generators\n')` and initialises a test result object via `new_test_result('kernel/pulses_waveform_suite', 'Pulse waveform generators', ...)`, where the stated requirement is that pulse waveform utilities must reproduce their analytic formulae and tabulated periodic sequences. Each check is performed with `test_close`, which appends explanatory messages to the result. The checks are:

- **Sawtooth and triangular waves**: with amplitude `amp=2`, frequency `freq=1`, and `time_grid=[0 0.25 0.5 0.75 1.0]`, `sawtooth` is compared against the reference `[-2 -1 0 1 -2]` (implementing `amplitude*(2*f*mod(t,1/f)-1)`), and `triwave` against the absolute value of that sawtooth reference. Tolerances are `1e-15` absolute and relative.
- **Uhrig delays**: for `T=1`, `N=3`, pulse positions follow `T*(sin(pi*(1:N)/(2*N+2)).^2-0.5)`; delays are adjacent differences, and the remaining time is split equally into a leading and trailing chunk. Compared with tolerances `1e-15`.
- **Periodic phase tables**: `pmlg5(1)` is compared against `(pi/180)*339.22` and `pmlg5(21)` against `pmlg5(1)` (20-pulse periodic table); `spinal(1)` against `(pi/180)*10` and `spinal(65)` against `spinal(1)` (64-pulse periodic table). Tolerances are `1e-14`.
- **Pulse envelopes**: `pulse_shape('rectangular',4)` against `ones(1,4)`; `pulse_shape('sinc3',3)` against `[0 pi 0]`, since sinc3 sampled at -3, 0, and +3 gives 0, pi, and 0. Tolerances are `1e-15`.
- **JCAMP pulse-file reading**: `read_wave('rectangular_1000.pk',4)` returns amplitude, phase, Cartesian controls, and a scaling factor; amplitude is compared against `ones(1,4)` (100 percent amplitude at every point), phase against `zeros(1,4)`, `Cx` against `ones(1,4)` (zero-phase polar coordinates convert into unit X control), `Cy` against `zeros(1,4)`, and the scaling factor against `1` (unit integral scaling). Tolerances are `1e-15`.
- **Veshtort-Griffin scaling**: `vg_pulse('E0A',7,2)` is compared against `vg_pulse('E0A',7,1)/2`, reflecting normalisation of VG pulse amplitudes as `2*pi*shape/duration`. Tolerances are `1e-14`.
- **WURST chirp**: `chirp_pulse(5,1,4,2,'wurst')` is evaluated on `time_grid=linspace(-0.5,0.5,5)`. References: amplitudes `2*pi*sqrt(4)*(1-abs(sin(pi*time_grid).^2))` (calibrated as `sqrt(bandwidth/duration)*(1-|sin(pi t)^p|)`), phases `pi*4*(time_grid.^2)` (linear chirp phase `pi*duration*bandwidth*t^2` on the normalised grid), frequencies `4*time_grid` (bandwidth times normalised time), durations `ones(1,5)/5` (uniform piecewise-constant slices), intervals `diff(time_grid)`, and Cartesian controls from `polar2cartesian(amps_ref,phis_ref)` (X is amplitude times cosine phase, Y is amplitude times sine phase). Amplitude, phase, frequency, and control comparisons use tolerances `1e-12`; duration and interval comparisons use `1e-15`.
- **Hyperbolic secant pulse**: `sech_pulse(3,2,5,2,5)` (peak amplitude 3, frequency modulation 2, phase modulation 5, duration 2, 5 points) is compared against a time grid `linspace(-dur/2,dur/2,npts)` centred at zero, amplitudes `peak_amp*sech(freq_mod*time_ref)`, phases `phase_mod*log(cosh(freq_mod*time_ref))`, and Cartesian controls from `polar2cartesian(amps_ref,phis_ref)`. Tolerances are `1e-15`.

## Inputs and outputs

### Syntax

```matlab
result=test_pulses_waveform_suite()
```

### Inputs

None. The function takes no arguments.

### Outputs

- `result` - regression test result object with explanatory messages, accumulated through repeated `test_close` calls.

## References

- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_pulses_waveform_suite.m>
- Spinach GitHub repository: <https://github.com/IlyaKuprov/Spinach>
