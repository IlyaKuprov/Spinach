# examples/nmr_spen/dosy_oneshot_1.m

- Signature: `dosy_oneshot_1()`

## Purpose

Oneshot DOSY pulse sequence for a system of three coupled spins with different relaxation rates. Timing: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Oneshot DOSY pulse sequence for a system of three coupled
- spins with different relaxation rates.
- Timing: minutes on NVidia Tesla A100, much longer on CPU
- Magnetic field
- Spin system
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Acquisition parameters
- Sample geometry
- Relaxation phantom
