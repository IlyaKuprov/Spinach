# experiments/imaging/spiral.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/imaging/spiral.m`
- Signature: `mri=spiral(spin_system,parameters,H,R,K,G,F)`
- Total lines: 211

## Purpose

2D imaging sequence with spiral sampling of the k-space. Syntax: mri=spiral_pulse_sequence(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.t_echo -echo time, seconds
- parameters.spiral_frq -angular frequency of the spiral, rad/s
- parameters.spiral_dur -duration of the spiral, seconds
- parameters.grad_amp -gradient amplitude at the end
- of the spiral, T/m

## Outputs

- mri -MRI image with square sinebell apodisation
- Prior to being Fourier transformed, the spiral is resampled onto
- a suitable square k-space grid -change the sequence to output
- spiral_kx, spiral_ky and spiral_z variables if you plan to pro-
- cess the data in a different way.

## Implementation structure

- 2D imaging sequence with spiral sampling of the k-space. Syntax:
- mri=spiral_pulse_sequence(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.t_echo -echo time, seconds
- parameters.spiral_frq -angular frequency of the spiral, rad/s
- parameters.spiral_dur -duration of the spiral, seconds
- parameters.grad_amp -gradient amplitude at the end
- of the spiral, T/m
- mri -MRI image with square sinebell apodisation
- Prior to being Fourier transformed, the spiral is resampled onto
- a suitable square k-space grid -change the sequence to output

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `step()`, `evolution()`, `spin()`, `cumsum()`, `spiral_x()`, `spiral_y()`, `kfigure()`, `ktitle()`, `kxlabel()`, `kylabel()`, `tic()`, `spiral_z()`, `toc()`.
