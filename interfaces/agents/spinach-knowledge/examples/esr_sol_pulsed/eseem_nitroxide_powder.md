# examples/esr_sol_pulsed/eseem_nitroxide_powder.m

- Signature: `eseem_nitroxide_powder()`

## Purpose

Powder-averaged two-pulse ESEEM on a 14N nitroxide radical. Time-domain simulation in Liouville space with powder averaging over a finite grid. Set to reproduce Figure 4a in http://dx.doi.org/10.1063/1.453532, ideal pulses are assumed. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Powder-averaged two-pulse ESEEM on a 14N nitroxide radical. Time-domain
- simulation in Liouville space with powder averaging over a finite grid.
- Set to reproduce Figure 4a in http://dx.doi.org/10.1063/1.453532, ideal
- pulses are assumed.
- Calculation time: seconds
- Magnet field
- System specification
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Set the sequence parameters
- Simulation
