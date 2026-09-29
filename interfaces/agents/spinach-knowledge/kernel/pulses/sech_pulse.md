# kernel/pulses/sech_pulse.m

- MATLAB source: [kernel/pulses/sech_pulse.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/sech_pulse.m)
- Spinach wiki: [sech_pulse.m](https://spindynamics.org/wiki/index.php?title=sech_pulse.m)
- Signature: `[Cx,Cy,time_grid,amps,phis]=sech_pulse(peak_amp,freq_mod,phase_mod,dur,npts)`

## Purpose

Constructs a sampled hyperbolic-secant RF pulse and returns Cartesian and amplitude-phase representations.

## Waveform

The source builds `time_grid=linspace(-dur/2,dur/2,npts)`, including the two duration endpoints when there is more than one point. At each point it computes `amps=peak_amp*sech(freq_mod*time_grid)` and `phis=phase_mod*log(cosh(freq_mod*time_grid))`, then calls `polar2cartesian(amps,phis)` to obtain `Cx` and `Cy`. The phase is zero at the pulse centre, `t=0`.

## Inputs and outputs

- `peak_amp` — real scalar peak amplitude, in radians per second.
- `freq_mod` — real scalar frequency-modulation parameter, in radians per second.
- `phase_mod` — real, dimensionless phase-modulation parameter.
- `dur` — positive real scalar duration, in seconds.
- `npts` — positive integer number of digitisation points.
- `Cx`, `Cy` — Cartesian RF coefficients, in radians per second.
- `time_grid` — sample times, in seconds; `amps` — sampled amplitudes, in radians per second; `phis` — sampled phases, in radians.

The source constructs these samples directly; it has no separate post-generation filter parameter.

## Example

```matlab
[Cx,Cy,time_grid]=sech_pulse(1,672,5,10.24e-3,1000);
plot(time_grid,[Cx; Cy]); kgrid; xlim tight;
kxlabel('time, seconds'); kylabel('amplitude, rad/s');
```
