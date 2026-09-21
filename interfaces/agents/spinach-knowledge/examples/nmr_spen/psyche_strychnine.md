# examples/nmr_spen/psyche_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/psyche_strychnine.m`
- Signature: `psyche_strychnine()`
- Total lines: 101

## Purpose

PSYCHE pure-shift NMR spectrum of strychnine. Calculation time: hours, faster on a GPU.

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- PSYCHE pure-shift NMR spectrum of strychnine.
- Calculation time: hours, faster on a GPU.
- Strychnine spin system
- Move shifts into spectral window
- Magnetic induction
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Saltire chirp parameters
- Coherent evolution timesteps
- Sample parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `state()`, `imaging()`, `fid()`, `fidps()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_2d()`, `plot_1d()`.
