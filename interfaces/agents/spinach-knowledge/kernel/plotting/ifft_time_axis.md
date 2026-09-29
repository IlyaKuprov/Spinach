# kernel/plotting/ifft_time_axis.m

- Signature: `[t_shift,t,dt]=ifft_time_axis(npts,df,zf)`

## Purpose

Builds unshifted and shifted time coordinates for an inverse FFT, with optional zero-fill points on both sides of the frequency-domain data. The function returns three outputs; `nifft` is an internal length, not a fourth output.

## Axis construction and units

The effective transform length is `nifft=npts+2*zf`, and `dt=1/(df*nifft)`. For `df` in Hz, `dt` and both time axes are in seconds.

- `t` is a column vector `(0:nifft-1)'*dt`, matching the unshifted `ifft` order.
- `t_shift` is the column vector `(-floor(nifft/2):ceil(nifft/2)-1)'*dt`, matching `fftshift(ifft(...))` order for either parity of `nifft`.
- Both axes have `nifft` samples; `dt` is a scalar.

## Inputs and guards

- `npts` - real numeric integer greater than 1.
- `df` - positive real numeric scalar frequency interval in Hz.
- `zf` - optional nonnegative real numeric integer; defaults to 0 and is added on each side, so the total length grows by `2*zf`.

Invalid inputs raise an error. There is no plotting or axes side effect; the routine only constructs coordinates.

## Links

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/ifft_time_axis.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ifft_time_axis.m)
