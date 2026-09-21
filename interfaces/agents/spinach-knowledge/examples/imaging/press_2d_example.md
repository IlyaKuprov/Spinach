# examples/imaging/press_2d_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/press_2d_example.m`
- Signature: `press_2d_example()`
- Total lines: 113

## Purpose

2D PRESS example. Three independent spin systems are localised in three spots of a 2D sample. The spots are slectively excited and their NMR spectra recorded. The followng are the frequencies to excite the three substances: parameters.rf_frq_list={-120e3 -100e3} -substance A parameters.rf_frq_list={-80e3 -10e3} -substance B parameters.rf_frq_list={+30e3 +100e3} -substance C Simulation time: minutes, faster with a Tes

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 2D PRESS example. Three independent spin systems are localised
- in three spots of a 2D sample. The spots are slectively excited
- and their NMR spectra recorded. The followng are the frequencies
- to excite the three substances:
- parameters.rf_frq_list={-120e3 -100e3} -substance A
- parameters.rf_frq_list={-80e3 -10e3} -substance B
- parameters.rf_frq_list={+30e3 +100e3} -substance C
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Magnetic induction
- Spin systems
- Basis set
- Disable path tracing

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `kfigure()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `ktitle()`, `imaging()`, `apodisation()`, `fftshift()`, `ifftshift()`, `plot_1d()`.
