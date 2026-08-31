# experiments/imaging/basic_1d_hard.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/basic_1d_hard.m`
- Signature: `fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F)`
- Total lines: 134

## Purpose

Basic 1D imaging sequence with a hard pulse. Syntax: fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.ro_grad_amp -readout gradient amplitude, T/m parameters.sweep -detection sweep width, Hz parameters.npoints -number of points in the fid parameters.offset -transmitter and receiver offset, Hz

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 31-32: Assemble the Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 34-35: Make pulse operators; implemented by `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 39-40: Hard 90-degree pulse; implemented by `parameters.rho0=step(spin_system,Hy,parameters.rho0,pi/2)`.
- Lines 42-44: Echo time; implemented by `parameters.rho0=evolution(spin_system,L,[],parameters.rho0, 0.5/parameters.sweep,parameters.npoints,'final')`.
- Lines 46-47: Hard 180-degree pulse; implemented by `parameters.rho0=step(spin_system,Hx,parameters.rho0,pi)`.
- Lines 49-51: Pre-phasing gradient; implemented by `parameters.rho0=evolution(spin_system,L-parameters.ro_grad_amp*G{1},[], parameters.rho0,0.5/parameters.sweep,parameters.npoints,'final')`.
- Lines 53-56: Acquisition under a gradient; implemented by `fid=acquire(spin_system,parameters,L+parameters.ro_grad_amp*G{1}, spalloc(size(L,1),size(L,2),0), spalloc(size(L,1),size(L,2),0))`.

### Key state/data transformations

- Lines 32: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 35: computes `Hp` using `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 36: computes `Hy` using `Hy=kron(speye(parameters.npts),(Hp-Hp')/2i)`.
- Lines 37: computes `Hx` using `Hx=kron(speye(parameters.npts),(Hp+Hp')/2)`.
- Lines 40: computes `parameters.rho0` using `parameters.rho0=step(spin_system,Hy,parameters.rho0,pi/2)`.
- Lines 54-56: computes `fid` using `fid=acquire(spin_system,parameters,L+parameters.ro_grad_amp*G{1}, spalloc(size(L,1),size(L,2),0), spalloc(size(L,1),size(L,2),0))`.

### Local helper functions

- Line 61: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Outputs

- fid -free induction decay that should be Fourier transformed
- to obtain the image

## Implementation structure

- Basic 1D imaging sequence with a hard pulse. Syntax:
- fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.ro_grad_amp -readout gradient amplitude, T/m
- parameters.sweep -detection sweep width, Hz
- parameters.npoints -number of points in the fid
- parameters.offset -transmitter and receiver offset, Hz
- fid -free induction decay that should be Fourier transformed
- to obtain the image
- Check consistency
- Assemble the Liouvillian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `evolution()`, `acquire()`, `spalloc()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `ischar()`, `isscalar()`.
