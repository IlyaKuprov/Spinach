# experiments/imaging/slice_select_1d.m

- Signature: `fid=slice_select_1d(spin_system,parameters,H,R,K,G,F)`

## Purpose

Slice selection diagnostics function. Executes a shaped pulse on the user-supplied 1D phantom and records a 1D image. Syntax: fid=slice_select_1d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.ss_grad_amp -the amplitude of slice selection
- gradient,T/m
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

- fid -k-space signal, run an FT to get the image

## Implementation structure

- Slice selection diagnostics function. Executes a shaped pulse on
- the user-supplied 1D phantom and records a 1D image. Syntax:
- fid=slice_select_1d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.ss_grad_amp -the amplitude of slice selection
- gradient,T/m
- parameters.rf_frq_list -a vector of RF frequencies at each
- pulse slice, Hz
- parameters.rf_amp_list -a vector of RF amplitudes at each
- pulse slice, rad/s
- parameters.rf_dur_list -a vector of pulse slice durations,
