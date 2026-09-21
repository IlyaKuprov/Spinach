# examples/nmr_spen/ufmq_4spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufmq_4spin.m`
- Signature: `ufmq_4spin()`
- Total lines: 114

## Purpose

4Q ultrafast MaxQ NMR spectrum for a coupled four-spin system in the presence of realistic diffusion. Calculation time: hours, much faster on GPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 4Q ultrafast MaxQ NMR spectrum for a coupled four-spin
- system in the presence of realistic diffusion.
- Calculation time: hours, much faster on GPU
- Magnetic field
- Chemical shifts
- 3J couplings
- 4J couplings
- 5J couplings
- Coherence selection
- Basis set
- Algorithmic options
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `spin()`, `imaging()`, `kfigure()`, `subplot()`, `fftshift()`, `plot_uf()`.
