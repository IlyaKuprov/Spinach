# examples/nmr_zerofield/earth_field_dfp.m

- Signature: `earth_field_dfp()`

## Purpose

Earth's field NMR Simulation for 2,6-difluoropyridine; replicates simulated spectra in Figure 7 of without the weighted addition of the uncoupled 1H signal. Calculation time: seconds.

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Earth's field NMR Simulation for 2,6-difluoropyridine; replicates
- simulated spectra in Figure 7 of
- without the weighted addition of the uncoupled 1H signal.
- Calculation time: seconds.
- Earth's field
- Isotopes and labels
- Chemical shifts
- J-couplings (experimental)
- Basis set
- Relaxation theory
- Sequence parameters
- This needs a GPU
