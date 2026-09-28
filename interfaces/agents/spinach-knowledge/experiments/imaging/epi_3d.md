# experiments/imaging/epi_3d.m

- Signature: `fid=epi_3d(spin_system,parameters,H,R,K,G,F)`

## Purpose

Diffusion-weighted 3D echo-planar imaging pulse sequence. Call it from the `imaging()` context, which provides `H`, `R`, `K`, `G`, and `F`.

## Parameters

- `parameters.image_size`: number of points in each dimension of the resulting image; the implementation requires a two-element integer vector with both values greater than one.
- `parameters.ss_grad_amp`: slice-selection gradient amplitude, T/m.
- `parameters.pe_grad_amp`: phase-encoding gradient amplitude, T/m.
- `parameters.pe_grad_dur`: phase-encoding gradient duration, s.
- `parameters.ro_grad_amp`: readout gradient amplitude, T/m.
- `parameters.ro_grad_dur`: readout gradient duration, s.
- `parameters.diff_g_amp` (optional): three-element vector of X, Y, and Z diffusion-gradient amplitudes, T/m, active during the echo-time intervals.
- `parameters.t_echo`: echo time after slice selection, s.
- `parameters.rho0`: initial state vector; `parameters.coil`: detection state vector.
- `parameters.npts`: three spatial grid sizes; `parameters.dims`: three sample dimensions.
- `parameters.rf_frq_list`, `parameters.rf_amp_list`, and `parameters.rf_dur_list`: equal-length shaped-pulse frequency, amplitude, and duration lists; `parameters.rf_phi`: pulse phase.

## Implementation

The background operator is `B=H+F+1i*R+1i*K`. A shaped proton pulse applies slice selection under `B+parameters.ss_grad_amp*G{1}`, followed by a rollback gradient. The state evolves for `parameters.t_echo`, receives an ideal 180-degree pulse, and evolves for another `parameters.t_echo`. If `parameters.diff_g_amp` is present, its three gradient components are included during both echo-time intervals.

After projecting the transverse proton signal for a three-dimensional spatial plot, the sequence applies phase-encoding and readout prephasers. It then precomputes propagators for phase encoding and alternating readout-gradient directions. Nested phase-encoding and readout loops detect with `parameters.coil` and fill `fid`. The loop runs on a GPU when `gpu` is enabled; the result is gathered before return.

## Output

- `fid`: k-space representation of the image.

## Reference

- [Spinach documentation for `epi_3d.m`](https://spindynamics.org/wiki/index.php?title=epi_3d.m)
