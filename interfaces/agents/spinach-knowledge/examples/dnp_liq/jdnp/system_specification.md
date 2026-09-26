# examples/dnp_liq/jdnp/system_specification.m

- Signature: `[sys,inter,bas,parameters]=system_specification()`

## Purpose

Parameters of the 2e1n system used for the simulations reported in https://doi.org/10.1039/d1cp04186j

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

## Implementation structure

- Parameters of the 2e1n system used for the simulations reported
- in https://doi.org/10.1039/d1cp04186j
- Spin system
- Nuclear chemical shift tensor
- Electron g-tensors -axial along Z
- Set coordinates
- Empty scalar coupling array
- Relaxation theory
- Basis set
- Tolerance settings
- Reference g-factors
