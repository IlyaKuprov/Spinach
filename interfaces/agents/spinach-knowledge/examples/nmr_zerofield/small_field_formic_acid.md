# examples/nmr_zerofield/small_field_formic_acid.m

- Signature: `small_field_formic_acid()`

## Purpose

Zero-field NMR spectroscopy -15N pyridine. Set to reproduce Fig 2 from http://dx.doi.org/10.1103/PhysRevLett.107.107601 Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Zero-field NMR spectroscopy -15N pyridine. Set to reproduce
- Fig 2 from http://dx.doi.org/10.1103/PhysRevLett.107.107601
- Calculation time: seconds
- Magnetic field
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
