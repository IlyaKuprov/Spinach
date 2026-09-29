# experiments/imaging/epi_3d.m

- Signature: `fid=epi_3d(spin_system,parameters,H,R,K,G,F)`.
- Canonical implementation: `experiments/imaging/epi_3d.m` — https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/epi_3d.m.

## Contract and physical sequence

Call from `imaging()`, which supplies `H`, `R`, `K`, `G`, and `F`. The function combines `B=H+F+1i*R+1i*K`. It applies a shaped RF pulse with slice selection along `G{1}` to the hard-coded `1H` channel, rolls back the slice gradient, then evolves for `t_echo`. Optional `diff_g_amp` adds a simultaneous three-axis diffusion-gradient vector during that interval; without it the interval is unweighted evolution. An ideal `Ly` 180-degree pulse then refocuses the transverse magnetisation.

The source projects the post-echo state into a voxel-resolved `1H` signal for display, then samples a two-dimensional EPI k-space grid: `G{2}` is stepped for phase encoding and `G{3}` is the alternating-polarity readout. The coil state detects each sample. The returned object is the sampled complex k-space array, not a reconstructed image. DNP polarisation is not prepared or measured here: `parameters.rho0` is caller-supplied, so no polarisation magnitude or image value is implied.

## Parameters and units

- `parameters.rho0`: column state vector with the same dimension as `H`.
- `parameters.coil`: detection state vector with the same dimension as `H`.
- `parameters.npts`: three positive integer voxel counts; `parameters.dims`: three finite positive spatial extents used for voxel display. This function does not state their unit.
- `parameters.image_size`: two odd finite integer counts, each at least 3, ordered as phase-encode rows and readout samples. The required `imaging()` context rejects even sizes before EPI runs.
- `parameters.ss_grad_amp`, `parameters.pe_grad_amp`, `parameters.ro_grad_amp`: real scalar gradient amplitudes in T/m for `G{1}`, `G{2}`, and `G{3}`, respectively.
- `parameters.pe_grad_dur`, `parameters.ro_grad_dur`: gradient durations in seconds; `parameters.t_echo`: positive finite echo interval in seconds.
- Optional `parameters.diff_g_amp`: three finite real amplitudes in T/m, applied along `G{1}`, `G{2}`, and `G{3}` during the echo interval.
- RF inputs `rf_frq_list`, `rf_amp_list`, and `rf_dur_list` are required with matching lengths; `rf_phi` is also required. Their units follow the RF pulse interface.
- `G` must contain at least three operators dimensionally compatible with `H`; `H`, `R`, `K`, and `F` must be same-size matrices.

## Returned data and source-derived numerical facts

`fid` is complex with dimensions `image_size(1)-by-image_size(2)`. The source samples readout with increments of `ro_grad_dur/(image_size(2)-1)` seconds and advances phase encoding with increments of `pe_grad_dur/(image_size(1)-1)` seconds. For example, the smallest accepted `image_size` is `3-by-3`, which allocates nine k-space samples and uses two intervals in each encoded direction; this is a shape example, not a simulated or measured result. The general increments above remain the source-defined timing formulas. The state is optionally moved to GPU when enabled and gathered before return.

## References

- [Spinach documentation: `epi_3d.m`](https://spindynamics.org/wiki/index.php?title=epi_3d.m).
- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/epi_3d.m).
