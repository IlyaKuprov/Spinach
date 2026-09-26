# experiments/nmr_liquids/pansy_cosy.m

- Signature: `fid=pansy_cosy(spin_system,parameters,H,R,K)`

## Purpose

Magnitude mode PANSY-COSY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=pansy_cosy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -nuclei on which the sequence runs,
- e.g. {'1H','13C'}
- parameters.sweep -a vector with two sweep widths in Hz
- parameters.npoints -a vector of integers specifying
- point count in each dimension
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.aa -magnitude mode COSY FID on F1,F1 nuclei
- fid.ab -magnitude mode COSY FID on F1,F2 nuclei
- Note: decoupling with respect to either of the working nuclei
- is impossible in this pulse sequence.
- Note: this is the magnitude-mode analytical-pathway version.
- Gradient echo/anti-echo PANSY-COSY variants should be
- implemented as separate pulse sequence functions.
- Andrew Porter

## Implementation structure

- Magnitude mode PANSY-COSY pulse sequence from:
- fid=pansy_cosy(spin_system,parameters,H,R,K)
- parameters.spins -nuclei on which the sequence runs,
- e.g. {'1H','13C'}
- parameters.sweep -a vector with two sweep widths in Hz
- parameters.npoints -a vector of integers specifying
- point count in each dimension
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.aa -magnitude mode COSY FID on F1,F1 nuclei
- fid.ab -magnitude mode COSY FID on F1,F2 nuclei
