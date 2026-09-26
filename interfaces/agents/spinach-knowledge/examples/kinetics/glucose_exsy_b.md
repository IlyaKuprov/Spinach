# examples/kinetics/glucose_exsy_b.m

- Signature: `glucose_exsy_b()`

## Purpose

2D EXSY of transmembrane exchange of 3,3-difluoroglucose. See the fitting example set for the script that yielded the para- meters used below. Calculation time: seconds

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 2D EXSY of transmembrane exchange of 3,3-difluoroglucose. See
- the fitting example set for the script that yielded the para-
- meters used below.
- Calculation time: seconds
- Magnet field
- Isotopes
- Chemical shifts
- J-couplings
- Cartesian coordinates
- Chemical subsystems
- Reaction rate matrix
- Equilibrium concentrations with alpha-beta imbalance
