# experiments/spen/idosyzs.m

- Signature: `inten=idosyzs(spin_system,parameters,H,R,K,G,F)`

## Purpose

A simplified model sequence of the ZS iDOSY pulse sequence. Syntax: inten=idosyzs(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H, R, K, G and F. Parameters: parameters.rho0 -initial state parameters.coil -detection state parameters.spins -nuclei on which the sequence runs parameters.g_amp -gradient amplitude for diffusion encoding (T/m) parameters.sel_

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Outputs

- inten -the absolute value of the first point in
- the free induction decay; this number is
- proportional to the integral of the real
- part of the correctly phased spectrum

## Implementation structure

- A simplified model sequence of the ZS iDOSY pulse sequence. Syntax:
- inten=idosyzs(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H, R, K, G and F. Parameters:
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.spins -nuclei on which the sequence runs
- parameters.g_amp -gradient amplitude for diffusion
- encoding (T/m)
- parameters.sel_g_amp -gradient amplitude during the
- selective pulse (T/m)
- parameters.rf_phi -phase of the inversion pulse (rad/s)
