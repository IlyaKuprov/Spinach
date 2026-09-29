# kernel/pulses/uhrig_times.m

MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/uhrig_times.m
Source Wiki page: https://spindynamics.org/wiki/index.php?title=uhrig_times.m

- Signature: `time_delays=uhrig_times(T,N)`

## Purpose

Returns the delay intervals for an Uhrig dynamical decoupling (UDD) sequence of `N` ideal pulses over total duration `T`.

## Timing construction

The implementation computes the centred pulse positions as `T*(sin(pi*(1:N)/(2*N+2)).^2-0.5)`, differences successive positions to obtain the interior intervals, and appends equal starting and trailing delays. Thus it returns `N+1` delay values: the first is before the first ideal pulse, the interior values are between pulses, and the last follows the final pulse. Their sum is `T`. The source comment attributes the timing formula to WSW's 2009 JCP paper.

## Parameters / inputs

- `T` - finite positive real scalar; total sequence duration, in seconds
- `N` - finite positive real integer scalar; number of pulses

## Output

- `time_delays` - `N+1` delays, in seconds
