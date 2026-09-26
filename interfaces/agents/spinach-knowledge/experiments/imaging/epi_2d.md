# experiments/imaging/epi_2d.m

- Signature: `mri=epi_2d(spin_system,parameters,H,R,K,G,F)`

## Purpose

Diffusion weighted echo planar 2D imaging pulse sequence with variable diffusion encoding direction. Syntax: mri=epi_2d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.pe_grad_dur -the duration of the phase encoding gradient (X), seconds parameters.ro_grad_dur -the duration of the readout gradient (Y), seconds parame

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Outputs

- mri -MRI image with square sinebell apodisation.

## Implementation structure

- Diffusion weighted echo planar 2D imaging pulse sequence with
- variable diffusion encoding direction. Syntax:
- mri=epi_2d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.pe_grad_dur -the duration of the phase encoding
- gradient (X), seconds
- parameters.ro_grad_dur -the duration of the readout
- gradient (Y), seconds
- parameters.image_size -number of points in each dimension
- of the resulting image
- parameters.diff_g_amp -[optional] a vector of diffusion
