# kernel/pulses/pulse_shape.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/pulse_shape.m
Spin Dynamics Wiki: https://spindynamics.org/wiki/index.php?title=pulse_shape.m

## Purpose

Generate one of four sampled real-valued pulse-envelope vectors. The function returns envelope samples only; it does not set pulse phase, absolute start time, duration, or physical amplitude units.

## Syntax

~~~matlab
waveform=pulse_shape(pulse_name,npoints)
~~~

- pulse_name: character vector naming gaussian, sinc5, sinc3, or rectangular.
- npoints: positive integer sample count.
- waveform: row vector of sampled envelope values.

## Implemented envelopes

- gaussian: evaluate normpdf(t)/sqrt(2) at npoints equally spaced samples of t from -2 to 2, inclusive.
- sinc5: evaluate pi*sinc(t) at equally spaced samples from -5 to 5, inclusive.
- sinc3: evaluate pi*sinc(t) at equally spaced samples from -3 to 3, inclusive.
- rectangular: return ones(1,npoints).

For the sinc cases, MATLAB's normalised sinc convention is used. The sampled coordinate is dimensionless in this function; no time vector or sampling interval is returned. The source formulas are used directly rather than rescaled to unit peak. The implementation's acceptance check is npoints>=1 with an integer requirement (although its error text says greater than 1).
