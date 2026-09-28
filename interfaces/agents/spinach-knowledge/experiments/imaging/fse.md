# experiments/imaging/fse.m

- Signature: `mri=fse(spin_system,parameters,H,R,K,G,F)`

## Purpose

Fast spin echo (FSE) imaging sequence. Call it from the `imaging()` context, which supplies `H`, `R`, `K`, `G`, and `F`.

## Physical / mathematical content

The sequence forms the background operator `B=H+F+1i*R+1i*K` and constructs a pulse operator for the first spin named in `parameters.spins`. It applies an initial 90-degree pulse, moves to the left edge of k-space under the readout gradient, then repeats a 180-degree pulse, phase-encoding gradient, readout, and opposite phase-encoding gradient for each image row.

## Numerical / algorithmic content

- Phase-encoding amplitudes span `-parameters.pe_grad_amp` to `parameters.pe_grad_amp` across `parameters.image_size(1)` rows.
- Each readout trajectory is generated with `krylov()` using a time step of `parameters.ro_grad_dur/(parameters.image_size(2)-1)`. Projection onto `parameters.coil` fills the corresponding k-space row.
- The k-space data receive square-sinebell apodisation in both dimensions. The output is `-real(fftshift(fft2(ifftshift(fid)),2))`.

## Parameters / inputs

- `spin_system`: must use the `sphten-liouv` or `zeeman-liouv` formalism.
- `H`, `R`, `K`, `F`: numeric matrices of the same dimensions.
- `G`: cell array containing at least two gradient operators; the sequence uses `G{1}` for phase encoding and `G{2}` for readout.
- `parameters.spins`: nonempty cell array of character strings; its first entry selects the spin used for pulses.
- `parameters.npts`: vector of positive integers used to construct the pulse operator.
- `parameters.rho0`: numeric initial state.
- `parameters.coil`: numeric detection operator.
- `parameters.pe_grad_amp`: real scalar phase-encoding gradient amplitude, T/m.
- `parameters.ro_grad_amp`: real scalar readout gradient amplitude, T/m.
- `parameters.pe_grad_dur`: positive real scalar phase-encoding gradient duration, seconds.
- `parameters.ro_grad_dur`: positive real scalar readout gradient duration, seconds.
- `parameters.image_size`: two integers greater than one specifying the number of points in each image dimension.

## Outputs

- `mri`: MRI image reconstructed from k-space data with square-sinebell apodisation.

## Reference

- <https://spindynamics.org/wiki/index.php?title=fse.m>