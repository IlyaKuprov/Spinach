# experiments/imaging/phase_enc_2d.m

- Signature: `mri=phase_enc_2d(spin_system,parameters,H,R,K,G,F)`.
- Canonical implementation: `experiments/imaging/phase_enc_2d.m` — https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/phase_enc_2d.m.

## Contract and physical sequence

Call from `imaging()`, which supplies `H`, `R`, `K`, `G`, and `F`. The background generator is `B=H+F+1i*R+1i*K`. The source applies a 90-degree pulse about `Ly` to `parameters.rho0`, evolves for `t_echo`, applies an ideal 180-degree pulse, then evolves for a second `t_echo`. If `diff_g_amp` is present, the same two-component gradient vector is applied during each echo interval; otherwise both intervals are free background evolution. Each phase-encode row is prepared with `G{1}`, pre-rolled under the reversed readout gradient `G{2}`, and sampled under the forward readout gradient with the coil state. Phase rows are independent and executed with `parfor`.

The source uses the first named spin in `parameters.spins`. The caller provides the starting state; this function neither prepares nor measures DNP polarisation.

## Parameters and units

- `parameters.rho0`: state vector matching `H`; `parameters.coil`: detection state vector matching `H`.
- `parameters.spins`: nonempty cell array of character strings; `parameters.npts`: vector of positive integer voxel counts.
- `parameters.image_size`: two odd integer counts, each at least three, ordered as phase rows and readout points; the `imaging()` context rejects even sizes.
- `parameters.t_echo`: positive echo interval in seconds. Optional `parameters.diff_g_amp` has two real gradient amplitudes in T/m along `G{1}` and `G{2}`.
- `parameters.pe_grad_amp`, `parameters.ro_grad_amp`: real scalar gradient amplitudes in T/m; `parameters.pe_grad_dur`, `parameters.ro_grad_dur`: corresponding durations in seconds.
- `G` must contain at least two operators. `H`, `R`, `K`, and `F` must be same-size matrices. The formalism must be `sphten-liouv` or `zeeman-liouv`.

## Returned image and source-derived numerical facts

The pre-reconstruction k-space array and returned image have dimensions `image_size(1)-by-image_size(2)`. The source spaces phase amplitudes with `linspace(-pe_grad_amp,+pe_grad_amp,image_size(1))` and samples readout at intervals of `ro_grad_dur/(image_size(2)-1)` seconds. It applies square-sinebell apodisation in both dimensions and returns the real 2D Fourier transform. For example, the minimum accepted `image_size` of `3-by-3` gives a 3-by-3 k-space array and output matrix; this is a dimension example only, not an image result. The source-defined phase spacing and readout timing are given above.

## References

- [Spinach documentation: `phase_enc_2d.m`](https://spindynamics.org/wiki/index.php?title=phase_enc_2d.m).
- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/phase_enc_2d.m).
