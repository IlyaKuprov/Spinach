# experiments/imaging/press_voxel_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/press_voxel_3d.m`
- Signature: `phan=press_voxel_3d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 193

## Purpose

Voxel selection diagnostics function for 3D PRESS sequences. Returns the sample excitation profile. Syntax: phan=press_voxel_3d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H, R, K, G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.ss_grad_amp -the three amplitudes of slice selection
- gradient, T/m
- parameters.rf_frq_list -cell array of three vectors of RF frequ-
- encies at each pulse slice, Hz
- parameters.rf_amp_list -cell array of three vectors of RF
- amplitudes at each pulse slice, rad/s
- parameters.rf_dur_list -cell array of three vectors of pulse
- slice durations, in seconds
- parameters.rf_phi -cell array of three pulse phases at
- time zero
- parameters.max_rank -cell array of three maximum rank in the
- Fokker-Planck pulse operator (2 is
- usually enough)

## Outputs

- phan -the excitation profile imprinted into a 3D phantom.
- Notes: add 'polyadic' to sys.enable, or this will crash your computer.

## Implementation structure

- Voxel selection diagnostics function for 3D PRESS sequences. Returns
- the sample excitation profile. Syntax:
- phan=press_voxel_3d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H, R, K, G, and F.
- parameters.ss_grad_amp -the three amplitudes of slice selection
- gradient, T/m
- parameters.rf_frq_list -cell array of three vectors of RF frequ-
- encies at each pulse slice, Hz
- parameters.rf_amp_list -cell array of three vectors of RF
- amplitudes at each pulse slice, rad/s
- parameters.rf_dur_list -cell array of three vectors of pulse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `polyadic()`, `opium()`, `state()`, `shaped_pulse_af()`, `evolution()`, `coherence()`, `fpl2phan()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `ischar()`, `any()`.
