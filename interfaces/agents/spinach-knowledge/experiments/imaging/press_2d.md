# experiments/imaging/press_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/press_2d.m`
- Signature: `fid=press_2d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 190

## Purpose

2D PRESS (voxel selective NMR) pulse sequence. Syntax: fid=press_2d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 47-48: Assemble the Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 50-51: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 55-58: Slice selection 90 degree pulse; implemented by `parameters.rho0=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(1)*G{1},Lx,Ly, parameters.rho0,parameters.rf_frq_list{1},parameters.rf_amp_list{1}, parameters.rf_du…`.
- Lines 60-62: Rephasing gradient; implemented by `parameters.rho0=evolution(spin_system,L-parameters.ss_grad_amp(1)*G{1},[], parameters.rho0,sum(parameters.rf_dur_list{1})/2,1,'final')`.
- Lines 64-66: Apply a crusher gradient; implemented by `parameters.rho0=evolution(spin_system,L+parameters.sp_grad_amp*(G{1}+G{2}),[], parameters.rho0,parameters.sp_grad_dur,1,'final')`.
- Lines 68-71: Slice selection 180 degree pulse; implemented by `parameters.rho0=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(2)*G{2},Lx,Ly, parameters.rho0,parameters.rf_frq_list{2},2*parameters.rf_amp_list{2}, parameters.rf_…`.
- Lines 77-78: Acquisition; implemented by `fid=acquire(spin_system,parameters,H+F,R,K)`.

### Key state/data transformations

- Lines 48: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 51: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 52: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 53: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 56-58: computes `parameters.rho0` using `parameters.rho0=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(1)*G{1},Lx,Ly, parameters.rho0,parameters.rf_frq_list{1},parameters.rf_amp_list{1}, parameters.rf_du…`.
- Lines 78: computes `fid` using `fid=acquire(spin_system,parameters,H+F,R,K)`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.ss_grad_amp -the two amplitudes of slice selection
- gradient, T/m
- parameters.rf_frq_list -cell array of two vectors of RF frequ-
- encies at each pulse slice, Hz
- parameters.rf_amp_list -cell array of two vectors of RF
- amplitudes at each pulse slice, rad/s
- parameters.rf_dur_list -cell array of two vectors of pulse
- slice durations, in seconds
- parameters.rf_phi -cell array of two pulse phases at
- time zero
- parameters.max_rank -cell array of two maximum rank in the
- Fokker-Planck pulse operator (2 is
- usually enough)
- parameters.sp_grad_amp -crusher gradient amplitude, T/m
- parameters.sp_grad_dur -crusher gradient duration, seconds

## Outputs

- fid -free induction decay of the NMR spectrum

## Implementation structure

- 2D PRESS (voxel selective NMR) pulse sequence. Syntax:
- fid=press_2d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.ss_grad_amp -the two amplitudes of slice selection
- gradient, T/m
- parameters.rf_frq_list -cell array of two vectors of RF frequ-
- encies at each pulse slice, Hz
- parameters.rf_amp_list -cell array of two vectors of RF
- amplitudes at each pulse slice, rad/s
- parameters.rf_dur_list -cell array of two vectors of pulse
- slice durations, in seconds

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `shaped_pulse_af()`, `evolution()`, `acquire()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `ischar()`, `any()`, `cellfun()`, `isvector()`, `isscalar()`.
