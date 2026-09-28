# kernel/optimcon/ctrl_trajan.m

- Signature: `ctrl_trajan(spin_system,waveform,traj_data,fidelities)`

## Purpose

Internal diagnostic plotting function for the optimal control module. It plots selected control-pulse and trajectory analyses. Specify plotting settings when setting up the optimal control problem in `optimcon.m`; do not call this function directly.

## Parameters / inputs

- `spin_system`: supplies plotting settings, pulse timing, integrator and control parameters, and the basis used for trajectory analysis.
- `waveform`: waveform supplied to a user-end function such as `grape_xy.m`; must be a real numeric array.
- `traj_data`: cell array of trajectory data structures returned by GRAPE, one per ensemble member.
- `fidelities`: fidelity array returned by a user-end function such as `grape_xy`; must be a real numeric array.

## Plotting and analysis

- Plot selections handled here include `spectrogram`, `xy_controls`, `phi_controls`, `amp_controls`, `frq_controls`, `correlation_order`, `coherence_order`, `local_each_spin`, `total_each_spin`, `level_populations`, and `robustness`. `time_by_slice` selects waveform slice number rather than elapsed seconds for the time axis. The `trajectory` key triggers a trajectory return; it does not produce a plot here.
- Cartesian control, phase, and amplitude plots use linear plots for the `trapezium` integrator and stairs plots for the `rectangle` integrator. Cartesian controls and paired-channel amplitudes are displayed in Hz; phases are displayed in radians. The Cartesian control plot shows lower and upper control bounds and reports the maximum nutation angle. The amplitude plot shows an amplitude bound.
- Spectrogram and instantaneous-frequency plots use paired control channels and require at least five initial slices of equal duration. Both use only that initial uniform-duration portion when later slice durations differ; the spectrogram labels its time axis as truncated in that case. Spectrogram, instantaneous-frequency, phase, and amplitude plots require an even number of control channels.
- Trajectory analyses process each ensemble member's forward trajectory with `trajan`, after tracing over spatial degrees of freedom. These analyses reject the `zeeman-hilb`, `zeeman-liouv`, and `zeeman-wavef` formalisms and request `sphten-liouv` instead.
- `robustness` plots a histogram of fidelities and displays their mean and standard deviation.

## Reference

- [Spinach documentation: ctrl_trajan.m](https://spindynamics.org/wiki/index.php?title=ctrl_trajan.m)