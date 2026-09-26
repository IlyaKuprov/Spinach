# examples/dnp_liq/odnp_liquid_2.m

- Signature: `odnp_liquid_2()`

## Purpose

Overhauser type DNP in liquid phase at room temperature, after a perfect inversion pulse on the electron ESR signal. The simulation uses Redfield theory to account for the dipolar cross-relaxation. Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Overhauser type DNP in liquid phase at room temperature, after a perfect
- inversion pulse on the electron ESR signal. The simulation uses Redfield
- theory to account for the dipolar cross-relaxation.
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Coordinates (Angstrom)
- Basis set
- Relaxation theory
- Spinach housekeeping
- Isotropic thermal equilibrium
- Electron control operator
