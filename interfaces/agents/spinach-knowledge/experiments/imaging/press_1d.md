# experiments/imaging/press_1d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/press_1d.m`
- Signature: `fid=press_1d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 152

## Purpose

1D PRESS (voxel selective NMR) pulse sequence. Syntax: fid=press_1d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 41-42: Assemble the Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 44-45: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 49-52: Slice selection pulse; implemented by `parameters.rho0=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp*G{1},Lx,Ly, parameters.rho0,parameters.rf_frq_list,parameters.rf_amp_list, parameters.rf_dur_list,pa…`.
- Lines 54-56: Rephasing gradient; implemented by `parameters.rho0=evolution(spin_system,L-parameters.ss_grad_amp*G{1},[], parameters.rho0,sum(parameters.rf_dur_list)/2,1,'final')`.
- Lines 58-59: Call to acquisition; implemented by `fid=acquire(spin_system,parameters,H+F,R,K)`.

### Key state/data transformations

- Lines 42: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 45: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 46: computes `Lx` using `Lx=kron(speye(prod(parameters.npts)),(Lp+Lp')/2)`.
- Lines 47: computes `Ly` using `Ly=kron(speye(prod(parameters.npts)),(Lp-Lp')/2i)`.
- Lines 50-52: computes `parameters.rho0` using `parameters.rho0=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp*G{1},Lx,Ly, parameters.rho0,parameters.rf_frq_list,parameters.rf_amp_list, parameters.rf_dur_list,pa…`.
- Lines 59: computes `fid` using `fid=acquire(spin_system,parameters,H+F,R,K)`.

### Local helper functions

- Line 64: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
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

- fid -free induction decay of the NMR spectrum

## Implementation structure

- 1D PRESS (voxel selective NMR) pulse sequence. Syntax:
- fid=press_1d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.ss_grad_amp -the amplitude of slice selection
- gradient, T/m
- parameters.rf_frq_list -a vector of RF frequencies at each
- pulse slice, Hz
- parameters.rf_amp_list -a vector of RF amplitudes at each
- pulse slice, rad/s
- parameters.rf_dur_list -a vector of pulse slice durations,
- in seconds

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `shaped_pulse_af()`, `evolution()`, `acquire()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `ischar()`, `isscalar()`, `isvector()`, `any()`.
