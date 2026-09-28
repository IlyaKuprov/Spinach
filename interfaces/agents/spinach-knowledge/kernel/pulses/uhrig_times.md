# kernel/pulses/uhrig_times.m

- Signature: `time_delays=uhrig_times(T,N)`

## Purpose

Returns the delay intervals for an Uhrig dynamical decoupling (UDD) sequence of `N` ideal pulses over total duration `T`.

## Numerical / algorithmic content

The interior pulse positions use the WSW 2009 JCP formula, `T*(sin(pi*(1:N)/(2*N+2)).^2-0.5)`. The function differences these positions to obtain the interior delays, then adds equal starting and trailing delays so the delay sum is `T`.

## Parameters / inputs

- `T` - finite positive real scalar; total sequence duration
- `N` - finite positive real integer; number of pulses

## Outputs

- `time_delays` - delays between ideal pulses; the first pulse follows the first delay, and a delay follows the last pulse

Source Wiki page: https://spindynamics.org/wiki/index.php?title=uhrig_times.m
