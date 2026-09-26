# kernel/plotting/ifft_time_axis.m

- Signature: `[t_shift,t,dt]=ifft_time_axis(npts,df,zf)`

## Purpose

Time axis for IFFT with optional zero-filling. Syntax: [t_shift,t,dt,nifft]=ifft_time_axis(npts,df,zf)

## Physical / mathematical content

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

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
