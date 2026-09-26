# examples/nmr_liquids/inept_strychnine.m

- Signature: `inept_strychnine()`

## Purpose

INEPT experiment on strychnine. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- INEPT experiment on strychnine.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Temperature
- Basis set
- Sequence parameters
- Create the spin system structure
- Generate isotopomers
- Preallocate the answer
- Loop over isotopomers
- Build the basis
