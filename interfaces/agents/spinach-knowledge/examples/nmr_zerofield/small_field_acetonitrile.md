# examples/nmr_zerofield/small_field_acetonitrile.m

- Signature: `small_field_acetonitrile()`

## Purpose

Small-field NMR spectroscopy -acetonitrile with 13C on the methyl group. Set to reproduce Figure 3 from Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Small-field NMR spectroscopy -acetonitrile with 13C on the
- methyl group. Set to reproduce Figure 3 from
- Calculation time: seconds
- Magnetic field, 2.64 mG
- Spin system
- Interactions
- Temperature
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
