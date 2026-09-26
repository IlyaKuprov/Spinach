# examples/nmr_zerofield/field_drop_acetonitrile.m

- Signature: `field_drop_acetonitrile()`

## Purpose

Zero-field NMR spectroscopy -acetonitrile. The simulation proceeds by computing the exact thermal equilibrium state and them propagating it through a time-dependent field drop. Set to reproduce Figure 7 from Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Zero-field NMR spectroscopy -acetonitrile. The simulation
- proceeds by computing the exact thermal equilibrium state
- and them propagating it through a time-dependent field drop.
- Set to reproduce Figure 7 from
- Calculation time: seconds
- Magnetic field (polariser)
- Spin system
- Interactions
- Temperature
- Basis set
- Sequence parameters
- Spinach housekeeping
