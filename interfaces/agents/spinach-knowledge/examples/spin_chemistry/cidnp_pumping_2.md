# examples/spin_chemistry/cidnp_pumping_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/cidnp_pumping_2.m`
- Signature: `cidnp_pumping_2()`
- Total lines: 84

## Purpose

A simulation of Figure 2A in IK's paper on chemically amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011). Calculation time: seconds.

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A simulation of Figure 2A in IK's paper on chemically amplified
- NOEs (https://doi.org/10.1016/j.jmr.2004.01.011).
- Calculation time: seconds.
- Magnet field
- Isotopes
- Chemical shifts
- Chemical shift anisotropies (DFT)
- Coordinates (DFT)
- J-coupling (expt)
- Relaxation theories
- Formalism and basis
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `relaxation()`, `state()`, `magpump()`, `unit_state()`, `evolution()`, `kfigure()`, `scale_figure()`, `subplot()`, `answer()`, `ktitle()`, `kxlabel()`.
