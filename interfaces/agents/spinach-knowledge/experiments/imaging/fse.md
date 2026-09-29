# experiments/imaging/fse.m

- Signature: `mri=fse(spin_system,parameters,H,R,K,G,F)`.
- Canonical implementation: `experiments/imaging/fse.m` — https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/fse.m.

## Contract and sequence

Call from `imaging()`, which supplies `H`, `R`, `K`, `G`, and `F`. The background generator is `B=H+F+1i*R+1i*K`. The first spin named in `parameters.spins` defines the transverse pulse operator. The code applies a 90-degree pulse about `Ly`, moves to the left edge of readout k-space with `G{2}`, and for each phase-encode row applies a 180-degree pulse, encodes with `G{1}`, samples the readout trajectory under `G{2}` through the coil state, then rewinds the phase gradient. In this implementation there is one 180-degree refocusing pulse and one readout per phase-encode row; it reconstructs after acquiring the rows.

The initial state `parameters.rho0` is caller-supplied. This is an MRI sequence, not a DNP preparation or polarisation measurement; no measured signal magnitude is reported.

## Parameters and units

- `parameters.rho0`: input state vector; `parameters.coil`: detection state vector. Both must match the dimension of `H`.
- `parameters.spins`: nonempty cell array of character strings; the first entry selects the observed spin.
- `parameters.npts`: vector of positive integer spatial sample counts; `parameters.image_size`: two odd integer counts, each at least three, [phase rows, readout points]; the `imaging()` context rejects even sizes.
- `parameters.pe_grad_amp`, `parameters.ro_grad_amp`: real scalar phase/readout gradient amplitudes in T/m along `G{1}` and `G{2}`.
- `parameters.pe_grad_dur`, `parameters.ro_grad_dur`: corresponding gradient durations in seconds.
- `G` must contain at least two operators. `H`, `R`, `K`, and `F` are same-size matrices. The source accepts `sphten-liouv` and `zeeman-liouv` formalisms.

## Returned image and source-derived numerical facts

The sampled complex array has dimensions `image_size(1)-by-image_size(2)`, with phase rows first and readout samples second. The source spaces phase amplitudes with `linspace(-pe_grad_amp,+pe_grad_amp,image_size(1))`; readout spacing is `ro_grad_dur/(image_size(2)-1)` seconds. It applies square-sinebell apodisation in both dimensions and returns a real matrix via `-real(fftshift(fft2(ifftshift(fid)),2))`. For example, the minimum accepted `image_size` of `3-by-3` produces a 3-by-3 sampled array and reconstructed matrix; this illustrates dimensions only, not an image result. The gradient spacing above is the source-defined timing formula.

## References

- [Spinach documentation: `fse.m`](https://spindynamics.org/wiki/index.php?title=fse.m).
- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/fse.m).
