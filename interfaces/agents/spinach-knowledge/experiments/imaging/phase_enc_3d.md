# experiments/imaging/phase_enc_3d.m

- Signature: `fid=phase_enc_3d(spin_system,parameters,H,R,K,G,F)`

## Purpose

3D MRI pulse sequence, with a slice selection stage followed by phase- encoded acquisition of the slice. Syntax: fid=phase_enc_3d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- parameters.ss_grad_amp -the amplitude of slice selection
- gradient,T/m
- parameters.pe_grad_amp -phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp -readout gradient amplitude, T/m
- parameters.ss_grad_dur -the duration of the slice selection
- gradient, seconds
- parameters.pe_grad_dur -the duration of the phase encoding
- gradient, seconds
- parameters.ro_grad_dur -the duration of the readout gradient,
- seconds
- parameters.image_size -number of points in each dimension of
- the resulting image

## Outputs

- fid -k-space representation of the image

## Implementation structure

- 3D MRI pulse sequence, with a slice selection stage followed by phase-
- encoded acquisition of the slice. Syntax:
- fid=phase_enc_3d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which would
- provide H,R,K,G, and F.
- parameters.ss_grad_amp - the amplitude of slice selection
- gradient,T/m
- parameters.pe_grad_amp - phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp - readout gradient amplitude, T/m
- parameters.ss_grad_dur - the duration of the slice selection
- gradient, seconds
- parameters.pe_grad_dur - the duration of the phase encoding
