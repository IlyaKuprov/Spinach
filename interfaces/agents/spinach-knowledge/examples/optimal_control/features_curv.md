# examples/optimal_control/features_curv.m

[Stable source link](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_curv.m)

## Calculation

The example transfers longitudinal magnetisation of two coupled 13C spins into their two-spin singlet state. At 14.1 T the chemical shifts are 0.00 and 0.25 ppm, and the scalar coupling is 60 Hz. The full `sphten-liouv` basis is used without approximation. The normalised initial state is `state(...,'Lz','all')`; the target is `singlet(spin_system,1,2)`. The drift is the NMR Hamiltonian, and the two RF controls are the 13C `Lx` and `Ly` operators on one channel.

The 50-slice waveform has 1 ms per slice. The power ensemble is the 11-point vector `2π×100×linspace(0.6,1.4,11)` as supplied in the source. A `SNS` penalty with weight 100 is enabled; the configured diagnostics are correlation order, coherence order, XY controls, robustness, and spectrogram. This is a simulated RF-power-robust control calculation; no experimental data are imported.

## Curvilinear GRAPE and output

The optimisation variables are amplitude and phase: the source maps `u=[r,φ]` to Cartesian controls `[r cos(φ), r sin(φ)]` and passes `dx_du` to `grape_curv`. The coded Jacobian matrix has rows `[cos(φ), sin(φ)]` and `[−r sin(φ), r cos(φ)]`. The initial guess uses unit amplitude and random phases scaled by π/2. `fmaxnewton` runs the curvilinear GRAPE objective; the optimised coordinates are converted to Cartesian controls and propagated by `shaped_pulse_xy` with `expv-pwc`. The final density operator is filtered to zero-quantum 13C coherence and two-spin correlation, then summarised by `stateinfo` (order 5).

The source header gives 1/√2 ≈ 0.7071 as the intended optimal fidelity, not as a value printed by this script. It describes ten power levels, while the executable `linspace(...,11)` setting supplies eleven. The results are simulations, not hardware measurements.
