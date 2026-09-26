# kernel/pulses/heterodyne.m

- Signature: `[X,Y]=heterodyne(dt,signal,freq)`

## Purpose

Signal heterodyne from wall clock time into the rotating frame using analytic signal demodulation: the negative frequency half of the spectrum (the counter-rotating component), as well as the direction-ambiguous DC and Nyquist bins, are dropped; the posi- tive frequency half is doubled and frequency-shifted. Syntax: [X,Y]=heterodyne(dt,signal,freq)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Parameters / inputs

- dt -time step in the input data, seconds
- signal -wall clock time signal, a column vector
- freq -frequency to be demodulated, Hz

## Outputs

- X, Y -in-phase and out-of-phase parts of the
- rotating frame signal, column vectors
- Notes: the signal must be sampled with more than two points per
- period of the frequency being demodulated; the transform
- is zero-phase, so sample k of the outputs refers to the
- same wall clock time as sample k of the input; dropping
- the DC bin subtracts the signal mean exactly; the record
- is treated as periodic by the FFT, and should therefore
- begin and end in dead time.

## Implementation structure

- Signal heterodyne from wall clock time into the rotating frame
- using analytic signal demodulation: the negative frequency half
- of the spectrum (the counter-rotating component), as well as the
- direction-ambiguous DC and Nyquist bins, are dropped; the posi-
- tive frequency half is doubled and frequency-shifted. Syntax:
- [X,Y]=heterodyne(dt,signal,freq)
- dt -time step in the input data, seconds
- signal -wall clock time signal, a column vector
- freq -frequency to be demodulated, Hz
- X, Y -in-phase and out-of-phase parts of the
- rotating frame signal, column vectors
- period of the frequency being demodulated; the transform
