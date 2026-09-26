# examples/kinetics/aziridine_exsy_1.m

- Signature: `aziridine_exsy_1()`

## Purpose

NOESY/EXSY experiment on phenylaziridine, including scalar relaxation of the second kind induced by the 14N nucleus, in a situation where the chemical exchange is relatively slow and scalar relaxation of the first kind does not manifest itself. Set to reproduce Figures 1b and 4b from All parameters, except for the isotropic chemical shifts, exchange ra- tes and correlation times come from a DFT calculation. Calculati

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- NOESY/EXSY experiment on phenylaziridine, including scalar relaxation
- of the second kind induced by the 14N nucleus, in a situation where
- the chemical exchange is relatively slow and scalar relaxation of the
- first kind does not manifest itself. Set to reproduce Figures 1b and
- 4b from
- All parameters, except for the isotropic chemical shifts, exchange ra-
- tes and correlation times come from a DFT calculation.
- Calculation time: minutes
- Magnet induction
- Isotopes
- Coordinates (Angstrom)
- 14N quadrupolar coupling
