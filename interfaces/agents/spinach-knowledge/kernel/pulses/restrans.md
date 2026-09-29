# kernel/pulses/restrans.m

- Signature: `[X,Y,dt]=restrans(X_user,Y_user,dt_user,omega,Q,model,up_factor)`

## Purpose

Applies the source's second-order RLC probe-response model to rotating-frame in-phase and out-of-phase waveform components and returns the filtered components on a coarser time grid.

## Time grid and signal processing

For `pwc`, each input pair is interpreted at a slice midpoint: the duration is `numel(X_user)*dt_user`, and the finer circuit grid is populated by nearest-neighbour interpolation with extrapolation. For `pwl` and `pwl_tsc`, the paired samples are interpreted at slice edges: the duration is `(numel(X_user)-1)*dt_user`, and linear interpolation fills the finer grid. The fine time step is `pi/(16*omega)` (the source labels this 16-fold Nyquist oversampling).

On that grid, the code forms the input carrier from amplitude `sqrt(X0.^2+Y0.^2)` and phase `atan2(Y0,X0)`, simulates the transfer function with numerator `1/Q` and denominator coefficients `[1/omega^2,1/(omega*Q),1]`, then heterodynes the response back into X and Y. Each component is demodulated with a `lowpass` call using an IIR impulse response, cutoff argument 1, sample-rate argument 64, and steepness 0.95. The samples are then downsampled using the integer stride obtained from `floor(numel(X)/(up_factor*numel(X_user)))`; `dt` is multiplied by that stride. For `pwl_tsc`, the internal time grid used by the diagnostic plot is shifted by `-2*Q/omega`; the function does not return a time grid.

## Inputs and outputs

`X_user` and `Y_user` must be real column vectors with equal lengths; `dt_user`, `omega`, and `Q` must be finite positive scalars. The model is one of `pwc`, `pwl`, or `pwl_tsc`, and `up_factor` is a finite positive integer. The code rejects `dt_user < pi/omega` as breaking its rotating-frame approximation. The source describes about 100 as a safe guess for `up_factor`; this is guidance in the source, not a general accuracy guarantee.

- `X`, `Y` — output rotating-frame components after the modeled response.
- `dt` — output slice duration in seconds.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/restrans.m) · [Spinach wiki page](https://spindynamics.org/wiki/index.php?title=restrans.m)
