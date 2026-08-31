# experiments/imaging/press_voxel_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/press_voxel_2d.m`
- Signature: `phan=press_voxel_2d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 173

## Purpose

Voxel selection diagnostics function for 2D PRESS sequences. Re- turns the sample excitation profile. Syntax: phan=press_voxel_2d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H, R, K, G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K,G,F)`.
- Lines 44-45: Compose Liouvillian; implemented by `L=H+F+1i*R+1i*K`.
- Lines 47-48: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 52-53: Override input with uniform initial condition; implemented by `Lz=state(spin_system,'Lz',parameters.spins{1})`.
- Lines 56-60: X slice selection pulse; implemented by `rho=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(1)*G{1},Lx,Ly, rho,parameters.rf_frq_list{1},parameters.rf_amp_list{1}, parameters.rf_dur_list{1},parameters.rf_…`.
- Lines 62-64: X rephasing gradient; implemented by `rho=evolution(spin_system,L-parameters.ss_grad_amp(1)*G{1},[], rho,sum(parameters.rf_dur_list{1})/2,1,'final')`.
- Lines 66-67: Isolate single-quantum; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},1}})`.
- Lines 69-73: Y slice selection pulse (scaled down to 90 degrees); implemented by `rho=shaped_pulse_af(spin_system,L+parameters.ss_grad_amp(2)*G{2},Lx,Ly, rho,parameters.rf_frq_list{2},parameters.rf_amp_list{2}, parameters.rf_dur_list{2}/2,parameters.r…`.
- Lines 75-77: Y rephasing gradient; implemented by `rho=evolution(spin_system,L-parameters.ss_grad_amp(2)*G{2},[], rho,sum(parameters.rf_dur_list{2})/4,1,'final')`.
- Lines 79-80: Isolate zero-quantum; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},0}})`.
- Lines 82-83: Get the phantom; implemented by `Lz=state(spin_system,'Lz',parameters.spins{1})`.

### Key state/data transformations

- Lines 45: computes `L` using `L=H+F+1i*R+1i*K`.
- Lines 48: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 49: computes `Lx` using `Lx=polyadic({{opium(prod(parameters.npts),1),(Lp+Lp')/2}})`.
- Lines 50: computes `Ly` using `Ly=polyadic({{opium(prod(parameters.npts),1),(Lp-Lp')/2i}})`.
- Lines 53: computes `Lz` using `Lz=state(spin_system,'Lz',parameters.spins{1})`.
- Lines 54: computes `rho` using `rho=kron(ones(prod(parameters.npts),1),Lz)`.
- Lines 84: computes `phan` using `phan=real(fpl2phan(rho,Lz,parameters.npts))`.

### Local helper functions

- Line 89: `grumble()` — `function grumble(spin_system,parameters,H,R,K,G,F)`.
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

## Outputs

- phan -the excitation profile imprinted into a 2D phantom.

## Implementation structure

- Voxel selection diagnostics function for 2D PRESS sequences. Re-
- turns the sample excitation profile. Syntax:
- phan=press_voxel_2d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H, R, K, G, and F.
- parameters.ss_grad_amp -the two amplitudes of slice selection
- gradient, T/m
- parameters.rf_frq_list -cell array of two vectors of RF frequ-
- encies at each pulse slice, Hz
- parameters.rf_amp_list -cell array of two vectors of RF
- amplitudes at each pulse slice, rad/s
- parameters.rf_dur_list -cell array of two vectors of pulse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `polyadic()`, `opium()`, `state()`, `shaped_pulse_af()`, `evolution()`, `coherence()`, `fpl2phan()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `ischar()`, `any()`.
