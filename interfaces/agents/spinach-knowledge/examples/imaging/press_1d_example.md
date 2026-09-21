# examples/imaging/press_1d_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/press_1d_example.m`
- Signature: `press_1d_example()`
- Total lines: 115

## Purpose

1D PRESS example. Three independent spin systems are localised in three areas of a 1D sample. The areas are slectively excited and their NMR spectra recorded. The followng are the frequencies to excite the three substances: pulse_frq=+100e3 -substances B and C pulse_frq=0; -substance C pulse_frq=-100e3 -substances A and C Simulation time: seconds, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1D PRESS example. Three independent spin systems are localised
- in three areas of a 1D sample. The areas are slectively excited
- and their NMR spectra recorded. The followng are the frequencies
- to excite the three substances:
- pulse_frq=+100e3 -substances B and C
- pulse_frq=0; -substance C
- pulse_frq=-100e3 -substances A and C
- Simulation time: seconds, faster with a Tesla V100 GPU.
- Magnetic induction
- Basis set
- Disable path tracing
- This needs a GPU

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `kfigure()`, `scale_figure()`, `subplot()`, `ylim()`, `ktitle()`, `kxlabel()`, `imaging()`, `apodisation()`, `fftshift()`, `ifftshift()`, `plot_1d()`.
