# examples/nmr_proteins/ct_hsqc_gb1.m

- Signature: `ct_hsqc_gb1()`
- Source: [examples/nmr_proteins/ct_hsqc_gb1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/ct_hsqc_gb1.m)

## Purpose

Simulates a two-dimensional constant-time HSQC experiment for GB1. The source estimates hours of simulation time and notes that a Tesla A100 GPU can make it faster.

## Physical / mathematical content

- Imports molecule 1 and backbone data from `2N9K.pdb` and `2N9K.bmrb`, with `noshift='delete'`; these are structural and shift inputs, not an experimental spectrum.
- Uses a 14.1 T field, interaction and proximity cutoffs of 2.0 and 4.0, and an IK-1 `sphten-liouv` basis with scalar-coupling connectivity, `inter_level=4`, and `prox_level=1`.
- Removes `13C` spins and sets `parameters.spins` to `15N` and `1H`, with `15N` decoupling in F2. This is a two-nucleus HSQC simulation, not a 3D triple-resonance sequence; the caller does not define a separate receiver-channel list.

## Numerical / algorithmic content

The caller sets `J=90`, sweeps [3000, 3000] Hz, offsets [-7300, 5100] Hz, 128 x 128 acquisition points, and 512 x 512 zero filling, with axes in ppm. The source does not state a unit for J. It simulates positive and negative FIDs with `liquid(...,@ct_hsqc,parameters,'nmr')`, applies squared-cosine apodisation to each, Fourier-transforms each along F2, and forms the States signal as the positive component plus the conjugated negative component. It then transforms F1 and plots the negative real part. No measured FID or spectrum is loaded.

## Implementation structure

- Sets the protein inputs, field, cutoffs, basis, greedy algorithm option, sequence parameters, and removes carbon spins.
- Builds the basis, simulates the two FID components, apodises and Fourier-transforms them, combines them, and plots `-real(spectrum)`.
