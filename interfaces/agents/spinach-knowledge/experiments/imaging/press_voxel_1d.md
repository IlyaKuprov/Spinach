# experiments/imaging/press_voxel_1d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/press_voxel_1d.m`
- Signature: `phan=press_voxel_1d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 154

## Purpose

Voxel selection diagnostics function for 1D PRESS sequences. Re- turns the sample excitation profile. Syntax: phan=press_voxel_1d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 42-43: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 45-46: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 50-51: Override input with uniform initial condition; implemented by `Lz=state(spin_system,'Lz',parameters.spins{1})`.
- Lines 54-58: Slice selection pulse; implemented by `rho=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp*G{1},Lx,Ly, rho,parameters.rf_frq_list,parameters.rf_amp_list, parameters.rf_dur_list,parameters.rf_phi, paramet…`.
- Lines 60-62: Rephasing gradient; implemented by `rho=evolution(spin_system,L-parameters.ss_grad_amp*G{1},[], rho,sum(parameters.rf_dur_list)/2,1,'final')`.
- Lines 64-65: Get the phantom; implemented by `phan=fpl2phan(rho,Lz,[parameters.npts 1])`.

### Key state/data transformations

- Lines 43: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 46: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 47: computes `Lx` using `Lx=polyadic({{opium(prod(parameters.npts),1),(Lp+Lp')/2}})`.
- Lines 48: computes `Ly` using `Ly=polyadic({{opium(prod(parameters.npts),1),(Lp-Lp')/2i}})`.
- Lines 51: computes `Lz` using `Lz=state(spin_system,'Lz',parameters.spins{1})`.
- Lines 52: computes `rho` using `rho=kron(ones(prod(parameters.npts),1),Lz)`.
- Lines 65: computes `phan` using `phan=fpl2phan(rho,Lz,[parameters.npts 1])`.

### Local helper functions

- Line 70: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.ss_grad_amp -the amplitude of slice selection
- gradient, T/m
- parameters.rf_frq_list -a vector of RF frequencies at each
- pulse slice, Hz
- parameters.rf_amp_list -a vector of RF amplitudes at each
- pulse slice, rad/s
- parameters.rf_dur_list -a vector of pulse slice durations,
- in seconds
- parameters.rf_phi -pulse phase at time zero
- parameters.max_rank -maximum rank in the Fokker-Planck
- pulse operator (2 is usually enough)

## Outputs

- phan -the excitation profile imprinted into a 1D phantom.

## Implementation structure

- Voxel selection diagnostics function for 1D PRESS sequences. Re-
- turns the sample excitation profile. Syntax:
- phan=press_voxel_1d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.ss_grad_amp -the amplitude of slice selection
- gradient, T/m
- parameters.rf_frq_list -a vector of RF frequencies at each
- pulse slice, Hz
- parameters.rf_amp_list -a vector of RF amplitudes at each
- pulse slice, rad/s
- parameters.rf_dur_list -a vector of pulse slice durations,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `polyadic()`, `opium()`, `state()`, `shaped_pulse_af()`, `evolution()`, `fpl2phan()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `ischar()`, `isscalar()`, `isvector()`.
