# experiments/imaging/epi_2d.m

- Signature: `mri=epi_2d(spin_system,parameters,H,R,K,G,F)`

## Purpose

Diffusion-weighted echo-planar 2D imaging pulse sequence with variable diffusion-encoding direction. Call this function from the `imaging()` context, which supplies `H`, `R`, `K`, `G`, and `F`.

## Parameters

- `parameters.pe_grad_dur`: phase-encoding gradient duration (X), in seconds; a positive real scalar.
- `parameters.ro_grad_dur`: readout gradient duration (Y), in seconds; a positive real scalar.
- `parameters.pe_grad_amp`, `parameters.ro_grad_amp`: phase-encoding and readout gradient amplitudes, each a real scalar.
- `parameters.image_size`: two integers greater than one, giving the number of points in each image dimension.
- `parameters.diff_g_amp` (optional): two real diffusion-gradient amplitudes in X and Y, in T/m. If specified, `parameters.diff_g_dur` is required.
- `parameters.diff_g_dur` (optional): diffusion-gradient duration in seconds; if present, a positive real scalar.
- `parameters.npts`: vector of positive integers used to construct the spatial pulse operators.
- `parameters.spins`: nonempty cell array of character strings; the first entry selects the spin for the pulse operators.
- `parameters.rho0`: numeric initial state.
- `parameters.coil`: numeric detection operator.

The function supports the `sphten-liouv` and `zeeman-liouv` formalisms. `H`, `R`, `K`, and `F` must be numeric matrices of the same size; `G` must be a cell array containing at least two gradient operators.

## Sequence and reconstruction

The background operator is `B=H+F+1i*R+1i*K`. The sequence applies an ideal 90-degree `Ly` pulse. If `diff_g_amp` is supplied, the X and Y diffusion gradients are applied for `diff_g_dur` before an ideal 180-degree `Lx` pulse and again after it.

Phase-encoding and readout gradients are prephased together for half the shorter gradient duration; any remaining half-duration of the longer prephaser is applied separately. The phase-encoding loop alternates the readout-gradient sign between lines, records the coil signal into k-space, and advances under the phase-encoding gradient. The k-space data receive square-sinebell apodisation in both dimensions.

## Output

- `mri`: MRI image computed as `real(fftshift(fft2(ifftshift(fid))))` from the apodised k-space data.

## Reference

<https://spindynamics.org/wiki/index.php?title=epi_2d.m>