# experiments/spen/ufmq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/ufmq.m`
- Signature: `fid=ufmq(spin_system,parameters,H,R,K,G,F)`
- Total lines: 262

## Purpose

Ultrafast multiple-quantum NMR, a literal implementation of Figure 1A from (http://dx.doi.org/10.1002/cphc.201800667). Syntax: fid=ufmq_nmr(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H, R, K, G, and F. Parameters: parameters.spins nuclei on which the sequence runs parameters.dims size of the sample, m parameters.npts number of grid points parameters.

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- fid -free induction decay of the ultrafast NMR spectrum.

## Implementation structure

- Ultrafast multiple-quantum NMR, a literal implementation of Figure 1A
- from (http://dx.doi.org/10.1002/cphc.201800667). Syntax:
- fid=ufmq_nmr(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H, R, K, G, and F. Parameters:
- parameters.spins nuclei on which the sequence runs
- parameters.dims size of the sample, m
- parameters.npts number of grid points
- parameters.npoints number of acquired points for each
- gradient readout
- parameters.nloops number of loop, where each loop consists of
- a positive and a negative readout

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `polyadic()`, `speye()`, `ismember()`, `inflate()`, `chirp_pulse()`, `step()`, `evolution()`, `coherence()`, `shaped_pulse_xy()`, `gpuArray()`, `report()`, `num2str()`, `fid()`, `gather()`.
