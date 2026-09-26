# experiments/nmr_liquids/ct_hsqc.m

- Signature: `fid=ct_hsqc(spin_system,parameters,H,R,K)`

## Purpose

Constant-time phase-sensitive HSQC pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Syntax

```matlab
fid=ct_hsqc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.decouple_f2 [optional] nuclei to decouple
- in F2, e.g. {'15N','13C'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.pos,fid.neg -two components of the States quadrature
- signal.
- Note: natural abundance simulations should make use of the isotope
- dilution functionality. See dilute.m function.

## Implementation structure

- Constant-time phase-sensitive HSQC pulse sequence from:
- fid=ct_hsqc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.decouple_f2 [optional] nuclei to decouple
- in F2, e.g. {'15N','13C'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.pos,fid.neg - two components of the States quadrature
