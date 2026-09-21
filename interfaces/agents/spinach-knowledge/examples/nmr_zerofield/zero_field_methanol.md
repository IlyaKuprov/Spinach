# examples/nmr_zerofield/zero_field_methanol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_zerofield/zero_field_methanol.m`
- Signature: `zero_field_methanol()`
- Total lines: 55

## Purpose

Zero-field NMR spectroscopy -13C methanol. Set to reproduce Figure 1 from http://dx.doi.org/10.1016/j.cplett.2013.06.042 Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Zero-field NMR spectroscopy -13C methanol. Set to reproduce
- Figure 1 from http://dx.doi.org/10.1016/j.cplett.2013.06.042
- Calculation time: seconds
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
