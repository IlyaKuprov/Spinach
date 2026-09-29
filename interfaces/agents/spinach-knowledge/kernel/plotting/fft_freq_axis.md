# kernel/plotting/fft_freq_axis.m

- Signature: `[f_shift,f,df]=fft_freq_axis(npts,dt,zf)`

Returns column-vector frequency axes for the unshifted and shifted ordering of an FFT, together with the bin spacing. No FFT is performed.

## Inputs

- `npts`: acquired-point count, an integer greater than 1.
- `dt`: positive real scalar time step.
- `zf`: non-negative integer number of added zero-fill points; defaults to zero when omitted.

## Construction

The transform length is `nfft=npts+zf`, the sampling frequency is `1/dt`, and `df=1/(dt*nfft). The unshifted axis is `f=(0:nfft-1)'*df`, with zero first. The shifted axis is `f_shift=(-floor(nfft/2):ceil(nfft/2)-1)'*df`, in the order corresponding to `fftshift(fft(...)). Both vectors have `nfft` rows; `nfft` is an internal length, not a returned output.

## Existing syntax

`[f_shift,f,df]=fft_freq_axis(npts,dt,zf)`

## Outputs

- `f_shift`: bins for data ordered as `fftshift(fft(...))`.
- `f`: bins for data ordered as `fft(...)`.
- `df`: frequency spacing.

## References

- [Source: `kernel/plotting/fft_freq_axis.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/fft_freq_axis.m)
- [Spinach Wiki: `fft_freq_axis.m`](https://spindynamics.org/wiki/index.php?title=fft_freq_axis.m)
