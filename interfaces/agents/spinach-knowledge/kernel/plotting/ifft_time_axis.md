# kernel/plotting/ifft_time_axis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/ifft_time_axis.m`
- Signature: `[t_shift,t,dt]=ifft_time_axis(npts,df,zf)`
- Total lines: 71

## Purpose

Time axis for IFFT with optional zero-filling. Syntax: [t_shift,t,dt,nifft]=ifft_time_axis(npts,df,zf)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Set default zero-fill; implemented by `if nargin<3, zf=0; end`.
- Lines 31-32: Check consistency; implemented by `grumble(npts,df,zf)`.
- Lines 34-35: IFFT length; implemented by `nifft=npts+2*zf`.
- Lines 37-38: Time resolution; implemented by `dt=1/(df*nifft)`.
- Lines 40-41: Build unshifted axis; implemented by `t=(0:(nifft-1)).'*dt`.
- Lines 43-44: Build shifted axis; implemented by `t_shift=(-floor(nifft/2):ceil(nifft/2)-1).'*dt`.

### Control flow inferred from the code

- Line 29: conditional branch on `nargin<3, zf=0; end`.

### Key state/data transformations

- Lines 35: computes `nifft` using `nifft=npts+2*zf`.
- Lines 38: computes `dt` using `dt=1/(df*nifft)`.
- Lines 41: computes `t` using `t=(0:(nifft-1)).'*dt`.
- Lines 44: computes `t_shift` using `t_shift=(-floor(nifft/2):ceil(nifft/2)-1).'*dt`.

### Local helper functions

- Line 49: `grumble()` — `function grumble(npts,df,zf)`.
  - Representative operation: `if (~isnumeric(npts))||(~isreal(npts))|| (~isscalar(npts))||(mod(npts,1)~=0)||(npts<=1)`.
  - Representative operation: `(~isscalar(npts))||(mod(npts,1)~=0)||(npts<=1)`.

## Parameters / inputs

- npts -number of frequency-domain points
- df -frequency interval between points, Hz
- zf -zero-fill length added to either
- side of the frequency domain

## Outputs

- t_shift -time axis for fftshift(ifft(...))
- t -time axis for ifft(...)
- dt -time step between points

## Implementation structure

- Time axis for IFFT with optional zero-filling. Syntax:
- [t_shift,t,dt,nifft]=ifft_time_axis(npts,df,zf)
- npts -number of frequency-domain points
- df -frequency interval between points, Hz
- zf -zero-fill length added to either
- side of the frequency domain
- t_shift -time axis for fftshift(ifft(...))
- t -time axis for ifft(...)
- dt -time step between points
- Set default zero-fill
- Check consistency
- IFFT length

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
