# examples/nmr_spen/ufmq_2spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufmq_2spin.m`
- Signature: `ufmq_2spin()`
- Total lines: 106

## Purpose

2Q ultrafast MaxQ NMR spectrum for a coupled two-spin system in the presence of realistic diffusion. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 2Q ultrafast MaxQ NMR spectrum for a coupled two-spin
- system in the presence of realistic diffusion.
- Calculation time: minutes on NVidia Tesla A100,
- much longer on CPU
- Magnetic field
- Chemical shifts
- J-coupling
- Coherence selection
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sample geometry

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `spin()`, `imaging()`, `kfigure()`, `subplot()`, `fftshift()`, `plot_uf()`.
