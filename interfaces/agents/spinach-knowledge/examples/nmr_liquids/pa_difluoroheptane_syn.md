# examples/nmr_liquids/pa_difluoroheptane_syn.m

- Signature: `pa_difluoroheptane_syn()`

## Purpose

Pulse-acquire 1H NMR spectrum of syn-3,5-difluoroheptane with a manual basis set specification as a merger of Lie algebras of the user-specified structral fragments followed by symmetry fac- torisation and conservation law screening. See our paper: for further information. Calculation time: minutes, faster with a GPU.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Pulse-acquire 1H NMR spectrum of syn-3,5-difluoroheptane with a
- manual basis set specification as a merger of Lie algebras of
- the user-specified structral fragments followed by symmetry fac-
- torisation and conservation law screening. See our paper:
- for further information.
- Calculation time: minutes, faster with a GPU.
- Magnet induction
- Isotopes
- Chemical shifts
- J-couplings
- Basis set
- GPU is useful here
