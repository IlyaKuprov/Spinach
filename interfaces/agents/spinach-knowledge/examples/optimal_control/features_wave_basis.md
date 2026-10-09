# examples/optimal_control/features_wave_basis.m

- Signature: `features_wave_basis()`

## Purpose

This example optimises coefficients of a pulse expressed in a user-selected Legendre basis. The source cites [doi:10.1016/j.jmr.2011.07.023](https://doi.org/10.1016/j.jmr.2011.07.023). Its design objective is excitation of the central 60 spins in a 100-spin offset array; the source explicitly leaves the dynamics of the 20 spins at either end unconstrained.

## Spin model and target

The model has 100 non-interacting `13C` spins at 14.1 T, with offsets equally spaced from −160 to +160 ppm. It uses `sphten-liouv` with `IK-2`, `prox_level=1` and `scalar_couplings` connectivity. The normalised optimisation states are `Lz` and `Lx` on spins 21–80, respectively.

## Basis-parameterised pulse design

The controls are `Lx` and `Ly` on ¹³C. The time grid has 125 slices of 4 μs (0.5 ms total); `wave_basis('legendre',20,125)'` supplies 20 basis coefficients per control channel. Eleven configured power levels span 15 × 10³ to 20 × 10³ rad/s after multiplication by 2π. The settings are `method='lbfgs'`, `max_iter=200`, and an `SNS` penalty of weight 10. The random initial coefficient array is `2 x 20`; the script calls `fmaxnewton` with `@grape_xy`, reconstructs the two-channel waveform from the coefficients and basis, then scales it by the mean configured power.

The later signal-generation section is distinct from the design objective: it initialises a new `Lz` state using the isotope label `13C`, simulates the shaped pulse, acquires with an `L+` coil, applies Gaussian apodisation, zero-fills and Fourier-transforms the FID, and plots the real spectrum. The acquisition settings are offset 0, sweep 55000, 2048 points and zero-fill 16384, with the axis labelled in ppm; the Gaussian apodisation argument is 10. The source does not annotate units for the sweep or apodisation value. The script requests a plot, but the source does not include the resulting spectrum or a convergence result, so neither outcome can be inferred from the setup.

Source: [examples/optimal_control/features_wave_basis.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_wave_basis.m).
