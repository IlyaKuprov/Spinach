# kernel/derivatives/fftdiff.m

- Signature: `kern=fftdiff(order,npoints,dx)`

## Purpose

Construct a Fourier-domain kernel for differentiating a periodic real signal sampled on a uniform grid. Apply it as:

```matlab
derivative=real(ifft(fft(signal).*kern));
```

The signal must have `npoints` samples with grid spacing `dx`. Periodic boundary conditions are assumed.

## Parameters / inputs

- `order` — positive integer derivative order.
- `npoints` — positive integer number of grid points.
- `dx` — positive real grid step length.

## Output

- `kern` — vector of Fourier-domain differentiation multipliers.

## Numerical / algorithmic content

Each kernel entry is the corresponding angular-frequency factor `(2i*pi*k/(npoints*dx))^order`. For odd `npoints`, the integer frequency indices run from `(1-npoints)/2` through `(npoints-1)/2`. For even `npoints`, they run from `-npoints/2` through `npoints/2-1`, including the negative Nyquist index. `ifftshift` places these multipliers in the ordering used by `fft`. The result is a spectral derivative of the periodic sampled signal, not a finite-difference approximation.

## Validation

The function rejects `order` or `npoints` unless each is a real numeric scalar integer of at least 1. It rejects `dx` unless it is a real numeric scalar greater than zero. The validation error messages describe `order` and `npoints` as “non-negative,” but the implemented checks require them to be positive.

## Reference

- [fftdiff.m](https://spindynamics.org/wiki/index.php?title=fftdiff.m)