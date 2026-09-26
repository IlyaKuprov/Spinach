# experiments/imaging/basic_1d_hard.m

- Signature: `fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F)`

## Purpose

Basic 1D imaging sequence with a hard pulse. Syntax: fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F. Parameters: parameters.ro_grad_amp -readout gradient amplitude, T/m parameters.sweep -detection sweep width, Hz parameters.npoints -number of points in the fid parameters.offset -transmitter and receiver offset, Hz

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Outputs

- fid -free induction decay that should be Fourier transformed
- to obtain the image

## Implementation structure

- Basic 1D imaging sequence with a hard pulse. Syntax:
- fid=basic_1d_hard(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F. Parameters:
- parameters.ro_grad_amp -readout gradient amplitude, T/m
- parameters.sweep -detection sweep width, Hz
- parameters.npoints -number of points in the fid
- parameters.offset -transmitter and receiver offset, Hz
- fid -free induction decay that should be Fourier transformed
- to obtain the image
- Check consistency
- Assemble the Liouvillian
