# kernel/pulses/oia_pulse.m

[Source: `kernel/pulses/oia_pulse.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/oia_pulse.m)

- Signature: `[Cx,Cy,durs,ints,amps,phis,frqs]=oia_pulse(npts,dur,bwidth,am_fun)`

## Purpose

Builds a frequency-swept inversion waveform with offset-independent adiabaticity (Tannus and Garwood, JMR A 120, 133 (1996)). The user supplies the amplitude modulation function; the frequency sweep is the normalised cumulative integral of the squared amplitude (Eq. 6 of the paper), so that the adiabaticity factor is the same at every offset inside the sweep bandwidth. The peak amplitude is calibrated so that the adiabaticity factor matches the inversion chirp of `chirp_pulse`; a constant amplitude function reproduces `chirp_pulse(npts,dur,bwidth,0,'smoothed')`.

## Inputs and discretisation

- `npts` — finite floating-point integer number of waveform points, at least 2.
- `dur` — finite positive floating-point pulse duration in seconds.
- `bwidth` — finite positive floating-point sweep bandwidth in Hz, centred on zero; integer-class scalars are refused.
- `am_fun` — function handle mapping a `1 x npts` row of normalised times in `[-1,1]` to a `1 x npts` row of non-negative finite floating-point amplitudes that are positive at every interior point (the envelope may vanish only at its two ends, because a zero inside the pulse also zeroes the sweep rate there); the scale does not matter because the envelope is normalised to unit peak internally.

The Table 1 shapes of the paper, with the 1% edge truncation used there, are `@(tau)1./(1+99*tau.^2)` (Lorentz), `@(tau)sech(asech(0.01)*tau)` (HS), `@(tau)exp(-log(100)*tau.^2)` (Gauss), `@(tau)(1+cos(pi*tau))/2` (Hanning), `@(tau)sech(asech(0.01)*tau.^n)` (HSn), and `@(tau)1-abs(sin(pi*tau/2)).^n` (Sin^n). Any smooth function that is positive inside the pulse and vanishes at most at its two ends is a valid envelope.

## Outputs and checks

- `Cx`, `Cy` — `1 x npts` real and imaginary RF components in rad/s.
- `durs` — `1 x npts` piecewise-constant slice durations in seconds; `ints` — `1 x (npts-1)` piecewise-linear interval durations in seconds.
- `amps`, `phis`, `frqs` — `1 x npts` amplitude in rad/s, phase in radians (zero at the centre of the pulse), and instantaneous frequency in Hz.

The routine validates the scalar inputs and the envelope samples, and refuses waveforms in which the phase advance across any interval exceeds π unless all seven outputs are requested.

## Reference

Tannus and Garwood, J. Magn. Reson. A 120, 133 (1996), Eq. 6 and Table 1. Regression: `tests/kernel/test_oia_pulse.m`; example: `examples/shaped_pulses/shaped_pulse_oia.m`.
