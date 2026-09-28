# experiments/imaging/phase_enc_2d.m

- Signature: `mri=phase_enc_2d(spin_system,parameters,H,R,K,G,F)`

## Purpose

2D phase-encoding imaging pulse sequence with optional diffusion weighting during the echo time. Call it from the `imaging()` context, which supplies `H`, `R`, `K`, `G`, and `F`.

## Sequence and reconstruction

The sequence forms the background operator `B=H+F+1i*R+1i*K` and applies an ideal 90-degree pulse about `Ly` to `parameters.spins{1}`. It evolves for `parameters.t_echo`, applies a 180-degree pulse about the same axis, and evolves for another `parameters.t_echo`. When `parameters.diff_g_amp` is present, both echo-time evolutions include `parameters.diff_g_amp(1)*G{1}+parameters.diff_g_amp(2)*G{2}`.

For each of `parameters.image_size(1)` phase-encoding steps, a gradient amplitude spanning `-parameters.pe_grad_amp` to `+parameters.pe_grad_amp` is applied with `G{1}` for `parameters.pe_grad_dur`. A pre-roll under `-parameters.ro_grad_amp*G{2}` precedes detection with `parameters.coil` under `+parameters.ro_grad_amp*G{2}`. The phase-encoding steps run in a `parfor` loop. The acquired data receive square-sinebell apodisation in both dimensions; the output is the real part of the shifted 2D Fourier transform.

## Parameters / inputs

- `spin_system`: spin system; its formalism must be `sphten-liouv` or `zeeman-liouv`.
- `H`, `R`, `K`, `F`: numeric matrices of the same dimensions, supplied by `imaging()`.
- `G`: cell array containing at least two gradient operators, supplied by `imaging()`.
- `parameters.npts`: vector of positive integers specifying spatial-grid point counts.
- `parameters.spins`: nonempty cell array of character strings; the first entry selects the spin for the pulse operator.
- `parameters.rho0`: numeric initial state.
- `parameters.coil`: numeric detection operator.
- `parameters.t_echo`: positive echo time in seconds, used on each side of the 180-degree pulse.
- `parameters.diff_g_amp`: optional two-element real vector of X- and Y-gradient amplitudes in T/m, active during the echo-time evolutions.
- `parameters.pe_grad_amp`: real phase-encoding gradient amplitude in T/m.
- `parameters.ro_grad_amp`: real readout gradient amplitude in T/m.
- `parameters.pe_grad_dur`: positive phase-encoding gradient duration in seconds.
- `parameters.ro_grad_dur`: positive readout gradient duration in seconds.
- `parameters.image_size`: two integers greater than one specifying the number of points in each image dimension.

## Output

- `mri`: real MRI image reconstructed with square-sinebell apodisation.

## Reference

- [Spinach documentation: phase_enc_2d.m](https://spindynamics.org/wiki/index.php?title=phase_enc_2d.m)