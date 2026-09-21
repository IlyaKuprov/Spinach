# experiments/imaging/fse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/fse.m`
- Signature: `mri=fse(spin_system,parameters,H,R,K,G,F)`
- Total lines: 178

## Purpose

Fast spin echo (FSE) pulse sequence. Syntax: mri=fse(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.pe_grad_amp -phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp -readout gradient amplitude, T/m
- parameters.pe_grad_dur -the duration of the phase encoding
- gradient, seconds
- parameters.ro_grad_dur -the duration of the readout gradient,
- seconds
- parameters.image_size -number of points in each dimension of
- the resulting image

## Outputs

- mri -MRI image with square sinebell apodisation.

## Implementation structure

- Fast spin echo (FSE) pulse sequence. Syntax:
- mri=fse(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.pe_grad_amp - phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp - readout gradient amplitude, T/m
- parameters.pe_grad_dur - the duration of the phase encoding
- gradient, seconds
- parameters.ro_grad_dur - the duration of the readout gradient,
- seconds
- parameters.image_size - number of points in each dimension of
- the resulting image

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `pe_grad_amps()`, `krylov()`, `rho_stack()`, `fid()`, `apodisation()`, `fftshift()`, `fft2()`, `ifftshift()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`.
