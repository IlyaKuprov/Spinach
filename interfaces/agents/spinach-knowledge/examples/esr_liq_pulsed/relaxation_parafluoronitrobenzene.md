# examples/esr_liq_pulsed/relaxation_parafluoronitrobenzene.m

- Signature: `relaxation_parafluoronitrobenzene()`

## Purpose

A pulse-acquire FFT version of the EasySpin parafluoronitrobenzene test file, with acknowledgements to Stefan Stoll. The Spinach simulation is run using explicit time propagation in Liouville space, including secular Redfield relaxation superoperator. Calculation time: minutes

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- A pulse-acquire FFT version of the EasySpin parafluoronitrobenzene test
- file, with acknowledgements to Stefan Stoll.
- The Spinach simulation is run using explicit time propagation in Liouville
- space, including secular Redfield relaxation superoperator.
- Calculation time: minutes
- Magnet field
- Isotope list
- Basis set
- Zeeman interactions
- Spin-spin couplings
- Relaxation superoperator
- Spinach housekeeping
