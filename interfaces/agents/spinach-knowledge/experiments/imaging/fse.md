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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 37-38: Assemble the background; implemented by `B=H+F+1i*R+1i*K`.
- Lines 40-41: Make pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 44-45: Apply 90-degree pulse; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 47-49: Move to left edge of k-space; implemented by `rho=step(spin_system,B+parameters.ro_grad_amp*G{2}, rho,parameters.ro_grad_dur/2)`.
- Lines 51-52: Preallocate k-space image; implemented by `fid=zeros(parameters.image_size,'like',1i)`.
- Lines 54-57: Phase encoding gradient amplitudes; implemented by `pe_grad_amps=linspace(-parameters.pe_grad_amp, parameters.pe_grad_amp, parameters.image_size(1))`.
- Lines 59-60: Get readout timing; implemented by `timestep=parameters.ro_grad_dur/(parameters.image_size(2)-1)`.
- Lines 62-63: Phase encoding loop; implemented by `for n=1:parameters.image_size(1)`.
- Lines 65-66: 180-degree pulse; implemented by `rho=step(spin_system,Ly,rho,pi)`.
- Lines 68-70: Kick up the k-space; implemented by `rho=step(spin_system,B+pe_grad_amps(n)*G{1}, rho,parameters.pe_grad_dur)`.
- Lines 72-75: Get trajectory using Krylov propagation; implemented by `rho_stack=krylov(spin_system,B+parameters.ro_grad_amp*G{2}, parameters.coil,rho,timestep, parameters.image_size(2)-1,'trajectory')`.
- Lines 78-80: Kick down the k-space; implemented by `rho=step(spin_system,B-pe_grad_amps(n)*G{1}, rho,parameters.pe_grad_dur)`.
- Lines 84-85: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}})`.
- Lines 87-88: Fourier transform; implemented by `mri=-real(fftshift(fft2(ifftshift(fid)),2))`.

### Control flow inferred from the code

- Line 63: `for` loop over `n=1:parameters.image_size(1)`.

### Key state/data transformations

- Lines 38: computes `B` using `B=H+F+1i*R+1i*K`.
- Lines 41: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 42: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 45: computes `rho` using `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 52: computes `fid` using `fid=zeros(parameters.image_size,'like',1i)`.
- Lines 55-57: computes `pe_grad_amps` using `pe_grad_amps=linspace(-parameters.pe_grad_amp, parameters.pe_grad_amp, parameters.image_size(1))`.
- Lines 60: computes `timestep` using `timestep=parameters.ro_grad_dur/(parameters.image_size(2)-1)`.
- Lines 73-75: computes `rho_stack` using `rho_stack=krylov(spin_system,B+parameters.ro_grad_amp*G{2}, parameters.coil,rho,timestep, parameters.image_size(2)-1,'trajectory')`.
- Lines 88: computes `mri` using `mri=-real(fftshift(fft2(ifftshift(fid)),2))`.

### Local helper functions

- Line 93: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

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
