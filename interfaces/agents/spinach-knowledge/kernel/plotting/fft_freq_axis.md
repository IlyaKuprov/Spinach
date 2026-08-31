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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Set default zero-fill; implemented by `if nargin<3, zf=0; end`.
- Lines 30-31: Check consistency; implemented by `grumble(npts,dt,zf)`.
- Lines 33-34: FFT length; implemented by `nfft=npts+zf`.
- Lines 36-37: Sampling frequency; implemented by `sfreq=1/dt`.
- Lines 39-40: Frequency resolution; implemented by `df=sfreq/nfft`.
- Lines 42-43: Build unshifted axis; implemented by `f=(0:(nfft-1)).'*df`.
- Lines 45-46: Build shifted axis; implemented by `f_shift=(-floor(nfft/2):ceil(nfft/2)-1).'*df`.

### Control flow inferred from the code

- Line 28: conditional branch on `nargin<3, zf=0; end`.

### Key state/data transformations

- Lines 34: computes `nfft` using `nfft=npts+zf`.
- Lines 37: computes `sfreq` using `sfreq=1/dt`.
- Lines 40: computes `df` using `df=sfreq/nfft`.
- Lines 43: computes `f` using `f=(0:(nfft-1)).'*df`.
- Lines 46: computes `f_shift` using `f_shift=(-floor(nfft/2):ceil(nfft/2)-1).'*df`.

### Local helper functions

- Line 51: `grumble()` — `function grumble(npts,dt,zf)`.
  - Representative operation: `if (~isnumeric(npts))||(~isreal(npts))|| (~isscalar(npts))||(mod(npts,1)~=0)||(npts<=1)`.
  - Representative operation: `(~isscalar(npts))||(mod(npts,1)~=0)||(npts<=1)`.

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
