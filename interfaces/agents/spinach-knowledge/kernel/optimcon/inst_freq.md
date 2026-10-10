# kernel/optimcon/inst_freq.m

- Signature: `freq=inst_freq(signal,dt,npoints,poly_order,amp_tol)`

## Purpose

Returns an instantaneous-frequency trajectory in hertz by unwrapping the phase of a complex signal and differentiating it with local Savitzky–Golay polynomial fits.

## Inputs and guards

All five arguments are required; the source assigns no defaults. `signal` must be a finite, non-empty complex numeric vector with at least three samples. Row and column vectors are accepted. `dt` must be a finite positive real scalar. `npoints` must be an odd integer of at least 3 and no greater than the sample count. `poly_order` must be a positive integer smaller than `npoints`. `amp_tol` must be a finite real scalar in [0,1].

## Calculation and output

The signal is columnised, its phase is formed with `unwrap(angle(signal))`, and `sgolaydiff(phase,1,npoints,poly_order)` is divided by `2*pi*dt`. At each sample, the routine positions a full `npoints`-sample stencil within the signal bounds, including at the edges. It sets the frequency to NaN if any sample in that stencil has magnitude less than or equal to `amp_tol*max(abs(signal))`. Thus `amp_tol=0` still masks stencils containing exact zero-magnitude samples. The output is reshaped to the input signal's row or column shape.

This routine has no optimiser, line search, freeze mask, or phase-cycle mask; its mask is solely the weak-amplitude stencil test above.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/inst_freq.m)
[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=inst_freq.m)
