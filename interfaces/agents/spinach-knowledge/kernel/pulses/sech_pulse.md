# kernel/pulses/sech_pulse.m

- Signature: `[Cx,Cy,time_grid,amps,phis]=sech_pulse(peak_amp,freq_mod,phase_mod,dur,npts)`

## Purpose

Constructs a hyperbolic-secant RF pulse and returns its Cartesian and amplitude-phase representations.

## Algorithm

The time grid is `linspace(-dur/2,dur/2,npts)`, centred on zero. The amplitude and phase are `peak_amp*sech(freq_mod*time_grid)` and `phase_mod*log(cosh(freq_mod*time_grid))`; `polar2cartesian` converts them to `Cx` and `Cy`.

## Parameters / inputs

- `peak_amp` — peak amplitude, in radians per second.
- `freq_mod` — frequency-modulation parameter, in radians per second.
- `phase_mod` — dimensionless phase-modulation parameter.
- `dur` — pulse duration, in seconds; must be positive.
- `npts` — number of digitisation points; must be a positive integer.

## Outputs

- `Cx`, `Cy` — coefficients of the `Sx` and `Sy` spin operators at each time slice, in radians per second.
- `time_grid` — pulse time points, in seconds.
- `amps` — pulse amplitudes, in radians per second.
- `phis` — pulse phases, in radians; the phase is zero at the centre (`t=0`).

## Example

```matlab
[Cx,Cy,time_grid]=sech_pulse(1,672,5,10.24e-3,1000);
plot(time_grid,[Cx; Cy]); kgrid; xlim tight;
kxlabel('time, seconds'); kylabel('amplitude, rad/s');
```

[Spinach wiki page](https://spindynamics.org/wiki/index.php?title=sech_pulse.m)
