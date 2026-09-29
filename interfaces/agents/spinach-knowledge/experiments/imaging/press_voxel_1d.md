# experiments/imaging/press_voxel_1d.m

Signature: `phan=press_voxel_1d(spin_system,parameters,H,R,K,G,F)`

## Contract and diagnostic profile

This imaging callback is normally invoked as imaging(spin_system,@press_voxel_1d,parameters). It is a diagnostic for the 1D PRESS selection profile, not an FID acquisition. The source accepts sphten-liouv or zeeman-liouv formalisms; H, R, K, and F must be same-sized matrices, and G must include G{1}.

The function creates a uniform longitudinal initial state by repeating Lz for spins{1} across npts sample points. It applies the shaped RF pulse with +ss_grad_amp G{1}, then evolves under the opposite gradient for half the sum of the RF durations. Finally it projects the state onto Lz with fpl2phan using grid shape [npts 1]. Thus the returned phan is a one-dimensional longitudinal/excitation profile with npts elements (an npts by 1 column for the scalar npts accepted here). It is not an acquired FID, a measured phantom, or an image reconstruction.

## Parameters and units

Required fields are spins (a nonempty cell array of spin labels), npts (positive integer scalar), ss_grad_amp (real scalar in T/m), rf_frq_list (Hz), rf_amp_list (rad/s), rf_dur_list (seconds), rf_phi (pulse phase), and positive-integer max_rank. The RF frequency, amplitude, and duration lists must have equal lengths. The function initialises its own uniform Lz state; it does not require the rho0 or coil fields used by the acquisition routine. max_rank controls the Fokker-Planck pulse operator; the source comment says 2 is usually enough.

## Source-backed configuration example

The shared example examples/imaging/press_1d_example.m uses npts 100 over a 0.30 m sample, ss_grad_amp 30 mT/m, rf_frq_list -100 kHz, rf_amp_list 2 pi times 5 kHz, rf_dur_list 50 microseconds, rf_phi pi/2, and max_rank 3. It supplies a three-component spatial phantom and longitudinal states for six 1H spins, then calls this function for the selection diagnostic. These are example inputs; no profile values are asserted here.

## References

- [Spinach Wiki: press_voxel_1d.m](https://spindynamics.org/wiki/index.php?title=press_voxel_1d.m)
- [Canonical MATLAB source: experiments/imaging/press_voxel_1d.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_voxel_1d.m)
- [Example configuration: press_1d_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/press_1d_example.m)
