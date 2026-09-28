# kernel/plotting/ifft_time_axis.m

- Signature: `[t_shift,t,dt]=ifft_time_axis(npts,df,zf)`

## Purpose

Constructs time axes for an inverse Fourier transform, accounting for optional zero padding on each side of the frequency domain. Syntax: `[t_shift,t,dt]=ifft_time_axis(npts,df,zf)`; the function returns three outputs (not `nifft`).

## Physical / mathematical content

- The time-step interval is `dt=1/(df*nifft)`, where `nifft=npts+2*zf`.

## Numerical / algorithmic content

- `t` is the column vector `(0:nifft-1).' * dt`; `t_shift` is `(-floor(nifft/2):ceil(nifft/2)-1).' * dt`. The function constructs these axes only; it does not modify data or perform an inverse Fourier transform.

## Parameters / inputs

- `npts` - number of frequency-domain points; real integer greater than 1
- `df` - frequency interval between points, in Hz; positive real scalar
- `zf` - zero-fill length on each side of the frequency domain; optional, defaults to 0, and must be a nonnegative real integer

## Outputs

- t_shift -time axis for fftshift(ifft(...))
- t -time axis for ifft(...)
- dt -time step between points