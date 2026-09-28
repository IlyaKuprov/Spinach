# kernel/plotting/fft_freq_axis.m

- Signature: `[f_shift,f,df]=fft_freq_axis(npts,dt,zf)`

## Purpose

Returns unshifted and FFT-shifted frequency axes and their frequency resolution for a time-domain acquisition with optional zero filling.

## Syntax

```matlab
[f_shift,f,df]=fft_freq_axis(npts,dt,zf)
```

## Parameters / inputs

- `npts` — number of acquired time-domain points; integer greater than 1.
- `dt` — time step between points; positive real scalar.
- `zf` — number of zero-fill points added to the acquired points; non-negative integer. Defaults to 0 when omitted.

## Numerical / algorithmic content

The transform length is `nfft=npts+zf`, the sampling frequency is `1/dt`, and the frequency resolution is `df=1/(dt*nfft)`. The unshifted axis `f` contains bins from 0 to `(nfft-1)*df`. The shifted axis is `(-floor(nfft/2):ceil(nfft/2)-1)*df`, matching the ordering produced by `fftshift`. The function validates the point count, time step, and zero-fill length.

## Outputs

- `f_shift` — frequency axis for data after `fftshift(fft(...))`.
- `f` — frequency axis for `fft(...)`.
- `df` — frequency resolution.
