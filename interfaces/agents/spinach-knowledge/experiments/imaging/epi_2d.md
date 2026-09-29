# experiments/imaging/epi_2d.m

- Signature: `mri=epi_2d(spin_system,parameters,H,R,K,G,F)`
- Canonical MATLAB source: [`experiments/imaging/epi_2d.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/epi_2d.m)

## Contract and sequence

This is a diffusion-weighted, two-dimensional echo-planar MRI simulation with variable diffusion-encoding direction. Call it from the `imaging()` context, which supplies `H`, `R`, `K`, `G`, and `F`. It forms `B=H+F+1i*R+1i*K`, uses `G{1}` along X for phase encoding and `G{2}` along Y for readout, and expands the requested `parameters.spins{1}` `Lx` and `Ly` pulse operators over the spatial grid.

The sequence starts with a hard pi/2 about `Ly`, optionally applies a diffusion gradient with X/Y amplitudes from `diff_g_amp`, applies a hard pi about `Lx`, and, when enabled, repeats that diffusion-gradient evolution. The prephasers bring the phase-encoding and readout gradients together for the shorter duration, with any remaining half-duration evolved separately. The readout propagators use opposite signs of the Y gradient. Successive phase-encoding rows alternate readout polarity; the negative-polarity row is stored in reverse readout order. After each row, the phase-encoding propagator advances the state.

The acquisition matrix `fid` has `parameters.image_size` shape. Each readout/phase increment is propagated for `ro_grad_dur/(image_size(2)-1)` or `pe_grad_dur/(image_size(1)-1)` seconds. The matrix is square-sinebell apodised in both dimensions and transformed as `real(fftshift(fft2(ifftshift(fid))))`; the returned `mri` therefore has the same two-dimensional grid size. These are simulated sampled data and a computed image, not a measured scan. The source does not specify a calibrated field of view or physical k-space units.

## Parameters, units, and constraints

- `parameters.pe_grad_dur` (phase encode, X) and `parameters.ro_grad_dur` (readout, Y): positive scalar durations in seconds.
- `parameters.image_size`: two odd integers, each at least 3; these set the phase and readout sample counts. The `imaging()` context rejects even sizes before EPI runs.
- `parameters.diff_g_amp`: optional two-element real X/Y vector in T/m; when present it requires `parameters.diff_g_dur`, a positive scalar in seconds.
- `parameters.pe_grad_amp` and `parameters.ro_grad_amp` are also used to scale `G{1}` and `G{2}` in the code. The source help text does not state their units, so use the gradient/operator convention of the calling imaging setup rather than assuming a calibration here.
- The function also requires numeric `parameters.rho0` and `parameters.coil`, a positive-integer `parameters.npts` vector, and nonempty `parameters.spins` cell array. `H`, `R`, `K`, and `F` must be same-size matrices, while `G` must contain at least two gradient operators.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=epi_2d.m>
- Source: [`experiments/imaging/epi_2d.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/epi_2d.m)
