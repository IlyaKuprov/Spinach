# examples/optimal_control/bloch_siegert/bloch_siegert_b.m

- Signature: `bloch_siegert_b()`
- Source: [examples/optimal_control/bloch_siegert/bloch_siegert_b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/bloch_siegert/bloch_siegert_b.m)

## Objective and spin model

Estimated calculation time: minutes. This example designs a universal rotation pulse over an ensemble of 100 non-interacting `13C` spins whose scalar chemical shifts are equally spaced from -100 to +100 ppm at `sys.magnet=28.18`. The selected basis is `sphten-liouv` with `IK-2`, proximity level 1, and scalar-coupling connectivity; the source comment says this retains each spin's complete basis while omitting multi-spin orders.

The target action is defined on three normalised states: `Sx -> -Sz`, `Sy -> Sy`, and `Sz -> Sx`. The control operators are `Lx` and `Ly` on `13C`, with the NMR-assumption Hamiltonian as drift. The three initial and target states are passed together to the ensemble optimiser.

## Waveform parameterisation and comparison

The 20 power levels span 0.001 to 1.0 times the absolute `13C` Larmor angular frequency and are specified in rad/s. Each waveform has 50 equal-duration slices, with `pulse_dt=(pi/100)/pwr_level`. A shared 2-by-50 Gaussian initial guess is scaled by 1/10. L-BFGS via `fmaxnewton` is limited to 500 iterations with `tol_x=1e-4`.

For every power, the code constructs optimiser settings with BSS off and on and optimises a `grape_xy` waveform under each. Both are then evaluated with the BSS-enabled settings using `ensemble`. The figure plots `1-fidelity` against relative power on a logarithmic scale.
