# examples/optimal_control/features_trapezium.m

- Signature: `features_trapezium()`

## Purpose

Optimise a pulse for state-to-state transfer across scalar couplings in a hydrofluorocarbon fragment spin system, from 1H Z-magnetisation to 19F Z-magnetisation. Six control operators define a piecewise-linear waveform using derivatives of a Lie group product quadrature, as described in [the cited paper](https://doi.org/10.1016/j.jmr.2023.107478). Calculation time: minutes.

## Physical / mathematical content

- The spin system contains 1H, 13C and 19F at a magnetic field of 9.4. All three chemical shifts are 0.0 ppm; scalar couplings are 140 Hz between spins 1 and 2 and -160 Hz between spins 2 and 3 (literature values).
- The basis uses `sphten-liouv` formalism with `none` approximation. The initial and target states are the normalised `Lz` states of spins 1 and 3, respectively.
- The drift operator is the NMR Hamiltonian. The six control operators are `Lx` and `Ly` for each of 1H, 13C and 19F.

## Numerical / algorithmic content

- The channel map is `[1; 1; 2; 2; 3; 3]`. Pulse-power levels are `2*pi*linspace(0.8e3,1.2e3,5)` rad/s; there are 50 slices of duration `2e-4`, with a waveform value at each of 51 endpoints.
- The control integrator is `trapezium`, the optimisation method is `lbfgs`, and `max_iter` is 100. The `SNS` penalty has weight 100.
- The initial guess is `randn(6,51)/3`. Optimisation calls `fmaxnewton(spin_system,@grape_xy,guess)`; the resulting waveform is multiplied by `mean(control.pwr_levels)`.
- A test simulation builds the left and right generators from successive waveform endpoints. Each slice is propagated with `step(spin_system,{G_L,(G_L+G_R)/2,G_R},rho,control.pulse_dt(n))`, retaining the piecewise-linear pulse model. The reported fidelity is `real(rho_targ'*rho)`.

## Implementation structure

- Create and basis-set the spin system, normalise the initial and target states, and construct the drift and control operators.
- Configure the controls and optimisation plots (`correlation_order`, `local_each_spin`, `xy_controls`, `spectrogram`), then call `optimcon` and optimise the pulse.
- Denormalise the pulse, run the test simulation and report `Re[<target|rho(T)>]`.
