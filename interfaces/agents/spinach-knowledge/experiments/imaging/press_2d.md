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
