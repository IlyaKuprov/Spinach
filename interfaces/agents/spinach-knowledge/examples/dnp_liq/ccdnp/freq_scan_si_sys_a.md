# examples/dnp_liq/ccdnp/freq_scan_si_sys_a.m

- Signature: `freq_scan_si_sys_a()`

## Purpose

Steady state nuclear magnetisation as a function of microwave frequency offset and the magnet field in a DNP experiment with two electrons con- nected by exchange coupling, both coupled to a nucleus by dipolar coup- lings. Further particulars in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Steady state nuclear magnetisation as a function of microwave frequency
- offset and the magnet field in a DNP experiment with two electrons con-
- nected by exchange coupling, both coupled to a nucleus by dipolar coup-
- lings. Further particulars in:
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Exchange coupling
- Coordinates for anisotropic HF
- Basis set
- Disable start-up checks
- Relaxation theory
