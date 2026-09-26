# examples/nmr_liquids/noesy_methanol.m

- Signature: `noesy_methanol()`

## Purpose

NOESY spectrum of 13C methanol. J-couplings from Pecul and Helgaker, CSA tensors from DFT. Note the presence of cross-peaks between 13C doublet components. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- NOESY spectrum of 13C methanol. J-couplings from Pecul and Helgaker,
- CSA tensors from DFT. Note the presence of cross-peaks between 13C
- doublet components.
- Calculation time: seconds
- Spin system properties (vacuum DFT calculation)
- Remove the OH proton
- Put all chemical shifts on resonance
- Magnet field
- Assign J-couplings
- Basis set
- Relaxation theory parameters
- Algorithmic options
