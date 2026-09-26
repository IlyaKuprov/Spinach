# experiments/imaging/udd_dec.m

- Signature: `mri=udd_dec(spin_system,parameters,H,R,K,G,F)`

## Purpose

The effect of Uhrid Dynamic Decoupling (UDD) pulse sequence on the MRI phantom. The function runs the UDD and then pro- jects out the user-specified spin state, returning the cor- responding image. Syntax: mri=udd_dec(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.dec_time -total duration of the sequence parameters.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Outputs

- mri -amplitude of the detection state at each point of the
- sample
- Note: the spin state to be observed should be specified in
- parameters.coil_st, the coil phantom is ignored.

## Implementation structure

- The effect of Uhrid Dynamic Decoupling (UDD) pulse sequence
- on the MRI phantom. The function runs the UDD and then pro-
- jects out the user-specified spin state, returning the cor-
- responding image. Syntax:
- mri=udd_dec(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.dec_time -total duration of the sequence
- parameters.npulses -number of pulses in the sequence,
- excluding the first pi/2 pulse
- parameters.spins -nuclei on which the sequence
- is to act, e.g. {'1H'}
