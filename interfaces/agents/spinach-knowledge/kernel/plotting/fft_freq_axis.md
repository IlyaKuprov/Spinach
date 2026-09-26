# kernel/plotting/fft_freq_axis.m

- Signature: `[f_shift,f,df]=fft_freq_axis(npts,dt,zf)`

## Purpose

Frequency axis for FFT with optional zero-filling. Syntax: [f_shift,f,df,nfft]=fft_freq_axis(npts,dt,zf)

## Physical / mathematical content

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

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
