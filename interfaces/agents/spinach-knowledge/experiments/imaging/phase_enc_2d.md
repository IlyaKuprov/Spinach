# experiments/imaging/phase_enc_2d.m

- Signature: `mri=phase_enc_2d(spin_system,parameters,H,R,K,G,F)`

## Purpose

2D phase encoding imaging pulse sequence with optional diffusion weighting during the echo time. Syntax: mri=phase_enc_2d(spin_system,parameters,H,R,K,G,F) This sequence must be called from the imaging() context, which would provide H,R,K,G, and F.

## Physical / mathematical content

- Imaging sequence implementations. They build spatially resolved Liouvillians that include gradients, slice-selection RF terms, diffusion, and acquisition operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- parameters.t_echo -echo time, seconds
- parameters.diff_g_amp -[optional] a vector of diffusion gra-
- dient pair amplitudes in X,Y (T/m) to
- be active during the echo time
- parameters.pe_grad_amp -phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp -readout gradient amplitude, T/m
- parameters.pe_grad_dur -the duration of the phase encoding
- gradient, seconds
- parameters.ro_grad_dur -the duration of the readout gradient,
- seconds
- parameters.image_size -number of points in each dimension of
- the resulting image

## Outputs

- mri -MRI image with square sinebell apodisation

## Implementation structure

- 2D phase encoding imaging pulse sequence with optional diffusion
- weighting during the echo time. Syntax:
- mri=phase_enc_2d(spin_system,parameters,H,R,K,G,F)
- This sequence must be called from the imaging() context, which
- would provide H,R,K,G, and F.
- parameters.t_echo - echo time, seconds
- parameters.diff_g_amp - [optional] a vector of diffusion gra-
- dient pair amplitudes in X,Y (T/m) to
- be active during the echo time
- parameters.pe_grad_amp - phase encoding gradient amplitude, T/m
- parameters.ro_grad_amp - readout gradient amplitude, T/m
- parameters.pe_grad_dur - the duration of the phase encoding
