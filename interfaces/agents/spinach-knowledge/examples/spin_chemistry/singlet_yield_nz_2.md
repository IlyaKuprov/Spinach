# examples/spin_chemistry/singlet_yield_nz_2.m

- Signature: `singlet_yield_nz_2()`

## Purpose

Field dependence of the decay rate of a micelle-confined triplet-born benzophenone ketyl / alkyl radical pair, computed with Redfield theory and with the lifetime-shifted Nakajima-Zwanzig kernel. The high-field decay of micellar pairs is relaxation-controlled: the T+/-states drain into the reactive S/T0 subspace at rates set by spectral densities of the anisotropic hyperfine modulation. Contact recombination at 1e9 H

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Field dependence of the decay rate of a micelle-confined triplet-born
- benzophenone ketyl / alkyl radical pair, computed with Redfield theory
- and with the lifetime-shifted Nakajima-Zwanzig kernel. The high-field
- decay of micellar pairs is relaxation-controlled: the T+/-states drain
- into the reactive S/T0 subspace at rates set by spectral densities of
- the anisotropic hyperfine modulation. Contact recombination at 1e9 Hz
- against a supercage correlation time of 0.7 ns puts the scalar lifetime
- shift at k*tau_c near 0.35, and the two theories separate. The observed
- decay rate is the slowest eigenmode of the pair Liouvillian. Parameters
- are representative of the SDS supercage systems of Sakaguchi, Hayashi,
- and Nagakura:
- Calculation time: minutes
