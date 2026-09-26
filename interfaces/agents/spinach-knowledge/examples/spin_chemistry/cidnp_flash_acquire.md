# examples/spin_chemistry/cidnp_flash_acquire.m

- Signature: `cidnp_flash_acquire()`

## Purpose

A model of the CIDNP magnetisation pumping process described in IK's paper: The system is pumped for 0.5 seconds, and then allowed to relax. Calculation time: seconds. Miguel Mompean Ilya Kuprov

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A model of the CIDNP magnetisation pumping process described in
- IK's paper:
- The system is pumped for 0.5 seconds, and then allowed to relax.
- Calculation time: seconds.
- Miguel Mompean
- Ilya Kuprov
- Magnet field
- Isotopes
- Chemical shifts
- Chemical shift anisotropies (DFT)
- Coordinates (DFT)
- J-coupling (expt)
