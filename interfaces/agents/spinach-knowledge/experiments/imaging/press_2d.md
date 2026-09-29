# experiments/imaging/press_2d.m

Signature: `fid=press_2d(spin_system,parameters,H,R,K,G,F)`

## Contract and selective-pulse sequence

This is an imaging callback, normally invoked as imaging(spin_system,@press_2d,parameters). The framework supplies H, R, K, F, gradient operators, initial state, and acquisition context. The source accepts sphten-liouv or zeeman-liouv formalisms; H, R, K, and F must be same-sized matrices, and the gradient cell must provide the two used directions G{1} and G{2}.

With L = H + F + iR + iK, the code constructs Lx and Ly from L+ for spins{1}, replicated over the spatial grid. First it applies a shaped 90-degree slice-selective RF pulse with ss_grad_amp(1) G{1}, rephases with the opposite G{1} gradient, and applies a crusher along G{1}+G{2}. It then applies the source-labelled 180-degree selective pulse along G{2}; the implementation doubles the second RF-amplitude list in this call. A second crusher is applied along G{1}+G{2}, after which the common acquire routine records the signal. The two selective directions define the active intersection; the crusher stages suppress unwanted coherence pathways in the simulated sequence.

## Parameters, units, and output

ss_grad_amp has two amplitudes, one per selection direction, in T/m. rf_frq_list, rf_amp_list, rf_dur_list, rf_phi, and max_rank are two-entry cell arrays, one per pulse; each pulse's RF frequency, amplitude, and duration lists must match in length. The source documents RF frequencies in Hz, amplitudes in rad/s, durations in seconds, and sp_grad_amp in T/m; sp_grad_dur is the crusher duration in seconds. spins is a nonempty cell array of spin labels, npts describes the spatial grid, and rho0 is the initial state. The downstream acquire routine also requires sweep (Hz), npoints, coil, and decouple. The returned fid is the one-dimensional complex acquisition FID with npoints samples, not a two-dimensional image. The example's image_size controls its associated voxel-selection diagnostic; it is not the PRESS FID shape.

## Source-backed configuration example

examples/imaging/press_2d_example.m sets ss_grad_amp to [25, 25] mT/m, crusher amplitude to 5 mT/m for 0.5 ms, pulse frequencies to -120 kHz and -100 kHz, amplitudes to 2 pi times 5 kHz, phases to pi/2, durations to 50 and 100 microseconds, and max_rank to 3 for each pulse. Its phantom grid is 108 by 90 over dimensions 0.30 by 0.25 m; image_size [101 105] is used for the separate 2D active-volume diagnostic. The example sets sweep 1000 Hz and npoints 128. These are configured values only, not an executed signal or measured result.

## References

- [Spinach Wiki: press_2d.m](https://spindynamics.org/wiki/index.php?title=press_2d.m)
- [Canonical MATLAB source: experiments/imaging/press_2d.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_2d.m)
- [Example configuration: press_2d_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/press_2d_example.m)
