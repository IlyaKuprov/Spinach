# kernel/pulses/heterodyne.m

- Signature: `[X,Y]=heterodyne(dt,signal,freq)`

## Purpose

Signal heterodyne from wall clock time into the rotating frame using analytic signal demodulation: the negative frequency half of the spectrum (the counter-rotating component), as well as the direction-ambiguous DC and Nyquist bins, are dropped; the positive frequency half is doubled and frequency-shifted.

## Parameters / inputs

- `dt` — time step in the input data, in seconds.
- `signal` — wall-clock-time signal, a column vector.
- `freq` — frequency to be demodulated, in hertz.

## Outputs

- `X`, `Y` — in-phase and out-of-phase components of the rotating-frame signal, as column vectors.

## Notes

- The signal must be sampled with more than two points per period of the frequency being demodulated.
- The transform is zero-phase, so output sample `k` refers to the same wall-clock time as input sample `k`.
- Dropping the DC bin subtracts the signal mean exactly.
- The FFT treats the record as periodic, so it should begin and end in dead time.
