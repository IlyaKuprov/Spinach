# experiments/imaging/phase_enc_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/phase_enc_2d.m`
- Signature: `mri=phase_enc_2d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 222

## Purpose

2D phase encoding imaging pulse sequence with optional diffusion weighting during the echo time. Syntax: mri=phase_enc_2d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 44-45: Assemble the background; implemented by `B=H+F+1i*R+1i*K`.
- Lines 47-48: Make pulse operators; implemented by `Hy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 51-52: Apply ideal 90-degree pulse; implemented by `parameters.rho0=step(spin_system,Hy,parameters.rho0,pi/2)`.
- Lines 54-55: Optional diffusion encoding; implemented by `if isfield(parameters,'diff_g_amp')`.
- Lines 57-60: Evolve for the echo time under a diffusion gradient; implemented by `parameters.rho0=step(spin_system,B+parameters.diff_g_amp(1)*G{1}+ parameters.diff_g_amp(2)*G{2}, parameters.rho0,parameters.t_echo)`.
- Lines 63-64: Just evolve for the echo time; implemented by `parameters.rho0=step(spin_system,B,parameters.rho0,parameters.t_echo)`.
- Lines 68-69: Apply 180-degree pulse; implemented by `parameters.rho0=step(spin_system,Hy,parameters.rho0,pi)`.
- Lines 85-88: Get PE gradient range; implemented by `pe_grad_amps=linspace(-parameters.pe_grad_amp, +parameters.pe_grad_amp, parameters.image_size(1))`.
- Lines 90-91: Preallocate the image; implemented by `fid=zeros(parameters.image_size,'like',1i)`.
- Lines 93-94: Phase encoding loop; implemented by `parfor n=1:parameters.image_size(1)`.
- Lines 96-97: Phase encoding gradient; implemented by `rho=evolution(spin_system,B+pe_grad_amps(n)*G{1},[],parameters.rho0,parameters.pe_grad_dur,1,'final')`.
- Lines 99-100: Get the timing parameters; implemented by `nsteps=parameters.image_size(2)-1; step_length=parameters.ro_grad_dur/nsteps`.
- Lines 102-103: Run the pre-roll gradient; implemented by `rho=evolution(spin_system,B-parameters.ro_grad_amp*G{2},[],rho,step_length/2,nsteps,'final')`.
- Lines 105-106: Detect under the readout gradient; implemented by `fid(n,:)=evolution(spin_system,B+parameters.ro_grad_amp*G{2},parameters.coil,rho,step_length,nsteps,'observable')`.
- Lines 110-111: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqsin'},{'sqsin'}})`.
- Lines 113-114: Fourier transform; implemented by `mri=real(fftshift(fft2(ifftshift(fid))))`.

### Control flow inferred from the code

- Line 55: conditional branch on `isfield(parameters,'diff_g_amp')`.
- Line 72: conditional branch on `isfield(parameters,'diff_g_amp')`.
- Line 94: `parfor` loop over `n=1:parameters.image_size(1)`.

### Key state/data transformations

- Lines 45: computes `B` using `B=H+F+1i*R+1i*K`.
- Lines 48: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 52: computes `parameters.rho0` using `parameters.rho0=step(spin_system,Hy,parameters.rho0,pi/2)`.
- Lines 86-88: computes `pe_grad_amps` using `pe_grad_amps=linspace(-parameters.pe_grad_amp, +parameters.pe_grad_amp, parameters.image_size(1))`.
- Lines 91: computes `fid` using `fid=zeros(parameters.image_size,'like',1i)`.
- Lines 97: computes `rho` using `rho=evolution(spin_system,B+pe_grad_amps(n)*G{1},[],parameters.rho0,parameters.pe_grad_dur,1,'final')`.
- Lines 100: computes `nsteps` using `nsteps=parameters.image_size(2)-1; step_length=parameters.ro_grad_dur/nsteps`.
- Lines 106: computes `fid(n,:)` using `fid(n,:)=evolution(spin_system,B+parameters.ro_grad_amp*G{2},parameters.coil,rho,step_length,nsteps,'observable')`.
- Lines 114: computes `mri` using `mri=real(fftshift(fft2(ifftshift(fid))))`.

### Local helper functions

- Line 119: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.t_echo -echo time, seconds
- parameters.diff_g_amp -[optional] a vector of diffusion gra-
- dient pair amplitudes in X,Y (T/m) to
- be active during the echo time
- parameters.pe_grad_amp -phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp -readout gradient amplitude, T/m
- parameters.pe_grad_dur -the duration of the phase encoding
- gradient, seconds
- parameters.ro_grad_dur -the duration of the readout gradient,
- seconds
- parameters.image_size -number of points in each dimension of
- the resulting image

## Outputs

- mri -MRI image with square sinebell apodisation

## Implementation structure

- 2D phase encoding imaging pulse sequence with optional diffusion
- weighting during the echo time. Syntax:
- mri=phase_enc_2d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.t_echo - echo time, seconds
- parameters.diff_g_amp - [optional] a vector of diffusion gra-
- dient pair amplitudes in X,Y (T/m) to
- be active during the echo time
- parameters.pe_grad_amp - phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp - readout gradient amplitude, T/m
- parameters.pe_grad_dur - the duration of the phase encoding

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `isfield()`, `evolution()`, `pe_grad_amps()`, `fid()`, `apodisation()`, `fftshift()`, `fft2()`, `ifftshift()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`.
