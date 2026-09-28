# examples/imaging/press_1d_example.m

- Signature: `press_1d_example()`

## Purpose

Demonstrates 1D PRESS localisation and spectral readout for three spin-pair components in a sample. The source describes the frequency settings used to scan the components and estimates seconds of runtime, faster with a Tesla V100 GPU.

## Spin systems and sample

At 3.0 T, six protons form three scalar-coupled pairs. Their chemical shifts are `{-3,+3}`, `{-2,+2}`, and `{-1,+1}` (the source's shift values), and the pair couplings are 10, 20, and 30. The basis uses the `sphten-liouv` formalism with `IK-2`, proximity level 1, and scalar-coupling connectivity. The 1D sample is 0.30 m long with 100 points; the initial-state phantoms define two localized regions and a third spatially uniform component.

## PRESS acquisition

The executed acquisition sets `rf_frq_list = -100e3`, RF phase `pi/2`, amplitude `2*pi*5000`, duration `0.5e-4 s`, and maximum rank 3. The source comments give the scan assignments as +100 kHz for B and C, 0 for C, and -100 kHz for A and C. A voxel-selection diagnostic is plotted before the PRESS signal is calculated. The signal receives square-cosine apodisation, then a magnitude Fourier transform is plotted as the voxel spectrum.
