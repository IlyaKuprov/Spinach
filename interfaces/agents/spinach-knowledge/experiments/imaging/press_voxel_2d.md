# experiments/imaging/press_voxel_2d.m

Source: [MATLAB on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_voxel_2d.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=press_voxel_2d.m)

- Signature: `phan=press_voxel_2d(spin_system,parameters,H,R,K,G,F)`

## Purpose and inputs

A parameterised 2D PRESS voxel-selection diagnostic, called from `imaging()` with the spatially resolved `H`, `R`, `K`, `G`, and `F` inputs. `parameters.ss_grad_amp` supplies two slice-gradient amplitudes in T/m. `parameters.rf_frq_list`, `parameters.rf_amp_list`, `parameters.rf_dur_list`, `parameters.rf_phi`, and `parameters.max_rank` are cell arrays with two entries, one per slice pulse: RF-frequency vectors in Hz, RF-amplitude vectors in rad/s, duration vectors in seconds, pulse phases at time zero (units are not specified by the source), and maximum Fokker–Planck pulse-operator ranks (2 is noted as usually sufficient), respectively. `parameters.spins{1}` identifies the selected spin; `parameters.npts` sets the spatial grid dimensions.

The source forms `L=H+F+1i*R+1i*K` and accepts `sphten-liouv` or `zeeman-liouv` formalism. It initialises a uniform `Lz` state across `prod(parameters.npts)` spatial points. This is a simulated diagnostic setup, not a measured initial profile.

## Sequence and returned data

The first shaped AFP pulse selects the X slice with `L+parameters.ss_grad_amp(1)*G{1}` and is rephased with the corresponding minus-gradient Liouvillian for half the sum of that pulse train's durations. The code then selects single-quantum coherence for `parameters.spins{1}`. The second shaped AFP pulse selects the Y slice at half the supplied second-pulse durations (the source labels this scaling as a 90-degree pulse); its rephasing evolution lasts one quarter of the second-pulse duration sum. The code selects zero-quantum coherence and calls `fpl2phan` with unweighted `coil_state` Lz as detection operator, distinct from the concentration-weighted initial `state`, returning `phan=real(...)` on the `parameters.npts` spatial grid. This function returns the 2D phantom array, not separate coordinate vectors.

This describes the implementation, not a run-verified or experimentally measured profile. No numerical phantom example or DOI is recorded in the source or existing page; the source's rank-2 note is retained above.
