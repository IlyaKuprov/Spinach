# kernel/pulses/heterodyne.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/heterodyne.m`
- Signature: `[X,Y]=heterodyne(dt,signal,freq)`
- Total lines: 83

## Purpose

Signal heterodyne from wall clock time into the rotating frame using analytic signal demodulation: the negative frequency half of the spectrum (the counter-rotating component), as well as the direction-ambiguous DC and Nyquist bins, are dropped; the posi- tive frequency half is doubled and frequency-shifted. Syntax: [X,Y]=heterodyne(dt,signal,freq)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Check consistency; implemented by `grumble(dt,signal,freq)`.
- Lines 40-41: Build time grid; implemented by `time_grid=dt*((1:numel(signal))'-1)`.
- Lines 43-44: One-sided spectral mask, DC and Nyquist bins dropped; implemented by `mask=zeros(numel(signal),1)`.
- Lines 47-48: Demodulate the analytic signal into the rotating frame; implemented by `signal=ifft(mask.*fft(signal)).*exp(-2i*pi*freq*time_grid)`.
- Lines 50-51: In-phase and out-of-phase components; implemented by `X=real(signal); Y=-imag(signal)`.

### Key state/data transformations

- Lines 41: computes `time_grid` using `time_grid=dt*((1:numel(signal))'-1)`.
- Lines 44: computes `mask` using `mask=zeros(numel(signal),1)`.
- Lines 45: computes `mask(2:ceil(numel(signal)/2))` using `mask(2:ceil(numel(signal)/2))=2`.
- Lines 48: computes `signal` using `signal=ifft(mask.*fft(signal)).*exp(-2i*pi*freq*time_grid)`.
- Lines 51: computes `X` using `X=real(signal); Y=-imag(signal)`.

### Local helper functions

- Line 56: `grumble()` — `function grumble(dt,signal,freq)`.
  - Representative operation: `if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))|| (~isfinite(dt))||(dt<=0)`.
  - Representative operation: `(~isfinite(dt))||(dt<=0)`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mask()`, `isscalar()`, `iscolumn()`.
