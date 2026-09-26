# examples/spin_chemistry/singlet_yield_nz_1.m

- Signature: `singlet_yield_nz_1()`

## Purpose

Magnetic field effect on a triplet-born benzophenone ketyl / thiyl radical pair in a viscous ionic liquid, computed with Redfield theory and with the lifetime-shifted Nakajima-Zwanzig kernel. The cage life- time of tens of nanoseconds and the rotational correlation time of a few nanoseconds put k*tau_c near 0.1, where the recombination drain visibly competes with rotational decorrelation of the anisotropic hyperfine 

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Magnetic field effect on a triplet-born benzophenone ketyl / thiyl
- radical pair in a viscous ionic liquid, computed with Redfield theory
- and with the lifetime-shifted Nakajima-Zwanzig kernel. The cage life-
- time of tens of nanoseconds and the rotational correlation time of a
- few nanoseconds put k*tau_c near 0.1, where the recombination drain
- visibly competes with rotational decorrelation of the anisotropic
- hyperfine couplings. Parameters are representative of the TMPA-TFSA
- measurements of the Wakasa group:
- Calculation time: minutes
- Field grid, cage drain, and correlation time
- Preallocate singlet yields
- Loop over the field grid
