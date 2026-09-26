# examples/nmr_liquids/pa_menthol.m

- Signature: `pa_menthol()`

## Purpose

Menthol NMR spectrum from Damien Jeannerat, including the effect of bad Z1 and Z2 magnet shims. Calculation time: minutes.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Menthol NMR spectrum from Damien Jeannerat, including the
- effect of bad Z1 and Z2 magnet shims.
- Calculation time: minutes.
- System and interaction specification
- Formalism and basis set
- Algorithms
- Spinach housekeeping
- Sequence parameters -1H
- Simulation
- Gaussian apodisation and then bad shims
- Fourier transform
- Plotting
