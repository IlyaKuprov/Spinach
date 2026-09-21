# kernel/plotting/fft_freq_axis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/fft_freq_axis.m`
- Signature: `[f_shift,f,df]=fft_freq_axis(npts,dt,zf)`
- Total lines: 74

## Purpose

Frequency axis for FFT with optional zero-filling. Syntax: [f_shift,f,df,nfft]=fft_freq_axis(npts,dt,zf)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- npts -number of acquired time-domain points
- dt -time step between points
- zf -zero-fill length added to time domain

## Outputs

- f_shift -frequency axis for fftshift(fft(...))
- f -frequency axis for fft(...)
- df -frequency resolution

## Implementation structure

- Frequency axis for FFT with optional zero-filling. Syntax:
- [f_shift,f,df,nfft]=fft_freq_axis(npts,dt,zf)
- npts -number of acquired time-domain points
- dt -time step between points
- zf -zero-fill length added to time domain
- f_shift -frequency axis for fftshift(fft(...))
- f -frequency axis for fft(...)
- df -frequency resolution
- Set default zero-fill
- Check consistency
- FFT length
- Sampling frequency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
