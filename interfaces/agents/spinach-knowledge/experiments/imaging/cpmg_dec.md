# experiments/imaging/cpmg_dec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/cpmg_dec.m`
- Signature: `mri=cpmg_dec(spin_system,parameters,H,R,K,G,F)`
- Total lines: 145

## Purpose

The effect of Carr-Purcell-Meiboom-Gill (CPMG) pulse sequence on the MRI phantom. The function runs the CPMG and then pro- jects out the user-specified spin state, returning the corres- ponding image. Syntax: mri=cpmg_dec(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.dec_time -total duration of the sequence paramet

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- mri -amplitude of the detection state at each point of the
- sample
- Note: the spin state to be observed should be specified in
- parameters.coil_st, the coil phantom is ignored.

## Implementation structure

- The effect of Carr-Purcell-Meiboom-Gill (CPMG) pulse sequence
- on the MRI phantom. The function runs the CPMG and then pro-
- jects out the user-specified spin state, returning the corres-
- ponding image. Syntax:
- mri=cpmg_dec(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.dec_time -total duration of the sequence
- parameters.npulses -number of pulses in the sequence,
- excluding the first pi/2 pulse
- parameters.spins -nuclei on which the sequence
- is to act, e.g. {'1H'}

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `report()`, `num2str()`, `fpl2phan()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `ischar()`, `isvector()`, `any()`, `isscalar()`.
