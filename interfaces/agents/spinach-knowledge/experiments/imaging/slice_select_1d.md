# experiments/imaging/slice_select_1d.m

Source: [MATLAB on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/slice_select_1d.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=slice_select_1d.m)

- Signature: `fid=slice_select_1d(spin_system,parameters,H,R,K,G,F)`

## Purpose and inputs

This parameterised sequence applies a shaped slice-selective pulse to the supplied 1D phantom state, rephases it, then delegates readout to `basic_1d_hard`. It is called from `imaging()` with `H`, `R`, `K`, `G`, and `F`; it requires `parameters.rho0` (numeric initial state), `parameters.ss_grad_amp` (slice-gradient amplitude in T/m), and shaped-pulse settings `parameters.rf_frq_list` (Hz), `parameters.rf_amp_list` (rad/s), `parameters.rf_dur_list` (seconds), `parameters.rf_phi` (pulse phase at time zero; units are not specified in this source), and `parameters.max_rank` (maximum Fokker–Planck pulse-operator rank; 2 is usually enough). `parameters.spins{1}` identifies the selected spin; `parameters.npts` defines the spatial grid over which the pulse acts. The called hard-readout helper also uses `parameters.ro_grad_amp` (readout-gradient amplitude, T/m), `parameters.sweep` (Hz), `parameters.npoints` (acquired FID point count), and `parameters.offset` (transmitter/receiver offset, Hz).

The sequence uses `L=H+F+1i*R+1i*K` and the source restricts it to `sphten-liouv` formalism. It constructs `Lx` and `Ly` over the `prod(parameters.npts)` spatial locations, applies the shaped AFP pulse with `L+parameters.ss_grad_amp*G{1}`, and rephases with the minus-gradient Liouvillian for half the sum of the RF-pulse durations. It acts on the caller-supplied `parameters.rho0`.

## Returned signal and axis interpretation

The function returns only `fid`, the k-space/free-induction signal produced by `basic_1d_hard`; Fourier transformation yields the image. The acquired signal length is controlled by `parameters.npoints`, while the phantom spatial discretisation is `parameters.npts`; this function does not return a separate coordinate-axis vector. This is a parameterised design, not a measured or run-verified image. No numeric worked example or DOI is present in the source or existing page; the source's rank-2 guidance is retained above.
