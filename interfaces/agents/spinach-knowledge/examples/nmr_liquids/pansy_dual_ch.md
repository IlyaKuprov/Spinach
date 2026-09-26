# examples/nmr_liquids/pansy_dual_ch.m

- Signature: `pansy_dual_ch()`

## Purpose

PANSY-COSY spectra of camphor with natural content of 13C isotope. Coordinates, shieldings, and J-couplings computed with DFT. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- PANSY-COSY spectra of camphor with natural content of 13C isotope.
- Coordinates, shieldings, and J-couplings computed with DFT.
- Calculation time: seconds
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Spinach housekeeping
- Sequence parameters
- Generate isotopomers
- Preallocate the answer
- Loop over isotopomers
- Build the basis
