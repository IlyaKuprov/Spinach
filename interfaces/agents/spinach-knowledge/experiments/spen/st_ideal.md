# experiments/spen/st_ideal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/st_ideal.m`
- Signature: `inten=st_ideal(spin_system,parameters,H,R,K,G,F)`
- Total lines: 132

## Purpose

The ideal Stejskal-Tanner pulse sequence using the notation from Figure 1 in http://dx.doi.org/0.1002/cmr.a.21241 with no gaps be- tween pulse sequence events. Syntax: inten=st_ideal(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H, R, K, G, and F. Parameters: parameters.spins -working spin. parameters.g_amp -gradient amplitude, T/m parameters.delta_sml 

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- inten -the absolute value of the first point in
- the free induction decay; this number is
- proportional to the integral of the real
- part of the correctly phased spectrum

## Implementation structure

- The ideal Stejskal-Tanner pulse sequence using the notation from
- Figure 1 in http://dx.doi.org/0.1002/cmr.a.21241 with no gaps be-
- tween pulse sequence events. Syntax:
- inten=st_ideal(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H, R, K, G, and F. Parameters:
- parameters.spins -working spin.
- parameters.g_amp -gradient amplitude, T/m
- parameters.delta_sml -the small delta parameter
- (see the figure)
- parameters.delta_big -the big delta parameter
- inten -the absolute value of the first point in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `elseif()`, `isfloat()`, `eps()`.
