# examples/nmr_spen/ufmq_6spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufmq_6spin.m`
- Signature: `ufmq_6spin()`
- Total lines: 125

## Purpose

6Q ultrafast MaxQ NMR spectrum for a coupled six-spin system in the presence of realistic diffusion. Calculation time: hours on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 6Q ultrafast MaxQ NMR spectrum for a coupled six-spin
- system in the presence of realistic diffusion.
- Calculation time: hours on NVidia Tesla A100, much longer on CPU
- Magnetic field
- Chemical shifts
- 3J couplings
- 4J couplings
- 5J couplings
- 6J couplings
- Coherence selection
- Basis set
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `spin()`, `imaging()`, `figure()`, `subplot()`, `fftshift()`, `plot_uf()`.
