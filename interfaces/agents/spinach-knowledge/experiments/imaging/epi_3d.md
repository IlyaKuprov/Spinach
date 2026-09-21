# experiments/imaging/epi_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/epi_3d.m`
- Signature: `fid=epi_3d(spin_system,parameters,H,R,K,G,F)`
- Total lines: 348

## Purpose

Diffusion weighted 3D echo planar imaging pulse sequence. Syntax: fid=epi_3d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.image_size -number of points in each dimension of the resulting image parameters.ss_grad_amp -the amplitude of slice selection gradient,T/m parameters.pe_grad_amp -phase encoding gradient ampli

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- fid -k-space representation of the image

## Implementation structure

- Diffusion weighted 3D echo planar imaging pulse sequence. Syntax:
- fid=epi_3d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.image_size - number of points in each dimension of
- the resulting image
- parameters.ss_grad_amp - the amplitude of slice selection
- gradient,T/m
- parameters.pe_grad_amp - phase encoding gradient amplitude, T/m
- parameters.pe_grad_dur - phase encoding gradient duration, s
- parameters.ro_grad_amp - readout gradient amplitude, T/m
- parameters.ro_grad_dur - readout gradient duration, s

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `shaped_pulse_af()`, `evolution()`, `isfield()`, `step()`, `fpl2phan()`, `state()`, `dims()`, `kfigure()`, `volplot()`, `ktitle()`, `propagator()`, `ismember()`, `gpuArray()`.
