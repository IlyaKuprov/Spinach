# experiments/imaging/epi_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/epi_2d.m`
- Signature: `mri=epi_2d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 214

## Purpose

Diffusion weighted echo planar 2D imaging pulse sequence with variable diffusion encoding direction. Syntax: mri=epi_2d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.pe_grad_dur -the duration of the phase encoding gradient (X), seconds parameters.ro_grad_dur -the duration of the readout gradient (Y), seconds parame

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 38-39: Assemble the background; implemented by `B=H+F+1i*R+1i*K`.
- Lines 41-42: Make pulse operators; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 47-48: Apply an ideal 90-degree pulse; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 50-51: Apply diffusion gradient; implemented by `if isfield(parameters,'diff_g_amp')`.
- Lines 57-58: Apply an ideal 180-degree pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 67-70: Preroll the gradients; implemented by `rho=evolution(spin_system,B-parameters.pe_grad_amp*G{1} -parameters.ro_grad_amp*G{2},[],rho, parameters.pe_grad_dur/2,1,'final')`.
- Lines 72-73: Preallocate k-space image; implemented by `fid=zeros(parameters.image_size,'like',1i)`.
- Lines 75-77: Precompute propagators; implemented by `P_ro_p=propagator(spin_system,B+parameters.ro_grad_amp*G{2}, parameters.ro_grad_dur/(parameters.image_size(2)-1))`.
- Lines 83-84: Phase encoding loop; implemented by `for n=1:parameters.image_size(1)`.
- Lines 86-87: Determine readout gradient sign; implemented by `ro_grad_sign=2*mod(n,2)-1`.
- Lines 89-90: Readout loop; implemented by `for k=1:parameters.image_size(2)`.
- Lines 92-93: Detect under readout gradient; implemented by `if ro_grad_sign>0`.
- Lines 103-104: Propagate under encoding gradient; implemented by `rho=P_pe_p*rho`.
- Lines 108-109: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}})`.
- Lines 111-112: Fourier transform; implemented by `mri=real(fftshift(fft2(ifftshift(fid))))`.

### Control flow inferred from the code

- Line 51: conditional branch on `isfield(parameters,'diff_g_amp')`.
- Line 61: conditional branch on `isfield(parameters,'diff_g_amp')`.
- Line 84: `for` loop over `n=1:parameters.image_size(1)`.
- Line 90: `for` loop over `k=1:parameters.image_size(2)`.
- Line 93: conditional branch on `ro_grad_sign>0`.

### Key state/data transformations

- Lines 39: computes `B` using `B=H+F+1i*R+1i*K`.
- Lines 42: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 43: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 48: computes `rho` using `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 73: computes `fid` using `fid=zeros(parameters.image_size,'like',1i)`.
- Lines 76-77: computes `P_ro_p` using `P_ro_p=propagator(spin_system,B+parameters.ro_grad_amp*G{2}, parameters.ro_grad_dur/(parameters.image_size(2)-1))`.
- Lines 78-79: computes `P_ro_m` using `P_ro_m=propagator(spin_system,B-parameters.ro_grad_amp*G{2}, parameters.ro_grad_dur/(parameters.image_size(2)-1))`.
- Lines 80-81: computes `P_pe_p` using `P_pe_p=propagator(spin_system,B+parameters.pe_grad_amp*G{1}, parameters.pe_grad_dur/(parameters.image_size(1)-1))`.
- Lines 87: computes `ro_grad_sign` using `ro_grad_sign=2*mod(n,2)-1`.
- Lines 94: computes `fid(n,k)` using `fid(n,k)=parameters.coil'*rho`.
- Lines 98: computes `fid(n,parameters.image_size(2)-k+1)` using `fid(n,parameters.image_size(2)-k+1)=parameters.coil'*rho`.
- Lines 112: computes `mri` using `mri=real(fftshift(fft2(ifftshift(fid))))`.

### Local helper functions

- Line 117: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Outputs

- mri -MRI image with square sinebell apodisation.

## Implementation structure

- Diffusion weighted echo planar 2D imaging pulse sequence with
- variable diffusion encoding direction. Syntax:
- mri=epi_2d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.pe_grad_dur -the duration of the phase encoding
- gradient (X), seconds
- parameters.ro_grad_dur -the duration of the readout
- gradient (Y), seconds
- parameters.image_size -number of points in each dimension
- of the resulting image
- parameters.diff_g_amp -[optional] a vector of diffusion

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `isfield()`, `evolution()`, `propagator()`, `fid()`, `apodisation()`, `fftshift()`, `fft2()`, `ifftshift()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`.
