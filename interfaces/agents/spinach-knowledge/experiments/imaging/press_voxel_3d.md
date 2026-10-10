# experiments/imaging/press_voxel_3d.m

Source: [MATLAB on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_voxel_3d.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=press_voxel_3d.m)

- Signature: `phan=press_voxel_3d(spin_system,parameters,H,R,K,G,F)`

## Purpose and inputs

A parameterised 3D PRESS voxel-selection diagnostic, called from `imaging()` with `H`, `R`, `K`, `G`, and `F`. `parameters.ss_grad_amp` contains three slice-gradient amplitudes in T/m. `parameters.rf_frq_list`, `parameters.rf_amp_list`, `parameters.rf_dur_list`, `parameters.rf_phi`, and `parameters.max_rank` are cell arrays with three entries, one per slice pulse: RF-frequency vectors in Hz, RF-amplitude vectors in rad/s, pulse-duration vectors in seconds, phases at time zero (units are not specified by the source), and maximum Fokker–Planck pulse-operator ranks (the source says 2 is usually enough). `parameters.spins{1}` identifies the selected spin, and `parameters.npts` defines the spatial grid.

The source forms `L=H+F+1i*R+1i*K`, accepts `sphten-liouv` and `zeeman-liouv`, and initialises a uniform `Lz` state across `prod(parameters.npts)` points. This is a simulated profile construction, not a measured voxel profile. The source comment says to add `polyadic` to `sys.enable`.

## Sequence and returned data

The X shaped AFP pulse uses the full first pulse train and is rephased for half its summed duration; the source then selects single-quantum coherence. The Y pulse uses half the second pulse durations (annotated as scaled to 90 degrees), is rephased for one quarter of their sum, and is followed by zero-quantum selection. The Z pulse likewise uses half the third pulse durations, is rephased for one quarter of their sum, and is followed by single-quantum selection. `fpl2phan` uses the unweighted `coil_state` L+ vector as detection operator, distinct from the concentration-weighted initial `state`; the returned `phan=abs(...)` is a 3D phantom array on `parameters.npts`, with no separate coordinate-vector outputs.
