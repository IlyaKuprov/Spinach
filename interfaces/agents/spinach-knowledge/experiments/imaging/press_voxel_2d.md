# experiments/imaging/press_voxel_2d.m

- Signature: `phan=press_voxel_2d(spin_system,parameters,H,R,K,G,F)`

## Purpose

Voxel selection diagnostics function for 2D PRESS sequences. Re- turns the sample excitation profile. Syntax: phan=press_voxel_2d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H, R, K, G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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
