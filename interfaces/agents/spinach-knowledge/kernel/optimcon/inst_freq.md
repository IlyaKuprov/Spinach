# kernel/optimcon/inst_freq.m

- Signature: `freq=inst_freq(signal,dt,npoints,poly_order,amp_tol)`

## Purpose

Computes the instantaneous frequency of a complex time-domain signal by unwrapping its phase and differentiating it with local Savitzky–Golay least-squares polynomial fits. The result is in hertz.

## Parameters / inputs

- `signal` — finite, non-empty complex vector with at least three samples.
- `dt` — positive sample interval in seconds.
- `npoints` — odd differentiation-window length, at least 3 and no greater than the signal length.
- `poly_order` — positive integer polynomial order smaller than `npoints`.
- `amp_tol` — amplitude threshold fraction from 0 to 1, relative to the signal's maximum magnitude.

## Output

- `freq` — instantaneous frequency in Hz, with the same shape and time grid as `signal`. A value is set to NaN if any sample in its local differentiation stencil has magnitude at or below the amplitude threshold; zero tolerance still masks zero-magnitude samples.

## Implementation

The routine unwraps `angle(signal)`, applies `sgolaydiff` to the phase, and divides by `2*pi*dt`. It then masks estimates whose local stencil contains a weak-signal sample.

[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=inst_freq.m)
