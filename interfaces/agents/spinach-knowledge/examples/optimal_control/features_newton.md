# examples/optimal_control/features_newton.m

- Signature: `features_newton()`

## Purpose

This example configures Newton-Raphson GRAPE pulse design for transfer from proton `Lz` to fluorine `Lz` across a three-spin H–C–F model. The stated design aim is robustness to proton transmitter offset and reduced pulse nutation frequency. The source cites Goodwin and Kuprov (2016), [doi:10.1063/1.4949534](https://doi.org/10.1063/1.4949534), for the method.

## Spin model and transfer

The model uses `1H`, `13C` and `19F` at 9.4 T, with all three chemical shifts set to 0.0 ppm. Its nonzero scalar couplings are 140 Hz for H–C and −160 Hz for C–F; the basis is `sphten-liouv` with approximation `none`. The normalised initial and target states are `Lz` on spin 1 (¹H) and spin 3 (¹⁹F), respectively.

## Configured pulse design

Six transverse controls (`Lx` and `Ly` for each isotope) share the corresponding three channels. The proton `Lz` operator defines five offset samples from −1000 to +1000 Hz. The configured pulse-power levels are `2*pi*[0.8 0.9 1.0]*1e3` rad/s. The waveform has 100 pointwise samples with 0.1 ms per slice (10 ms total); an `NS` penalty is assigned weight 0.01. The source describes a penalty when the waveform exceeds a user-specified power threshold, but does not give a separate numeric threshold.

The script sets `method='newton'` and `max_iter=50`, starts from a random `6 x 100` guess, and calls `fmaxnewton(...,@grape_xy,guess)`. These are design settings, not evidence that optimisation converged. Afterward it runs one shaped-pulse simulation under the configured drift and reports the real target overlap; the source contains no resulting fidelity value and does not show separate post-optimisation simulations over the offset/power ensemble.

Source: [examples/optimal_control/features_newton.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_newton.m).