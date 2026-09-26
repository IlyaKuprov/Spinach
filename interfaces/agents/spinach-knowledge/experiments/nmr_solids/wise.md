# experiments/nmr_solids/wise.m

- Signature: `fid=wise(spin_system,parameters,H,R,K)`

## Purpose

WISE (WIdeline SEparation) is a powder MAS heteronuclear correlation experiment. In the common 1H-13C implementation, molecular dynamics information is contained in 1H line shapes that are separated in the second dimension by 13C chemical shifts. Further information in:

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Syntax

```matlab
fid=wise(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -working spins, e.g. {'1H',13C'}
- parameters.hi_pwr -amplitude of high power pulses
- on the high-gamma channel, Hz
- parameters.cp_pwr -amplitude of CP pulse on each
- channel during the CP contact
- time, Hz
- parameters.cp_dur -CP contact time duration, s
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.sweep -sweep width, Hz for F1, F2
- parameters.npoints -number of points in F1, F2
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.sin, fid.cos -sine and cosine components
- of the States quadrature

## Implementation structure

- WISE (WIdeline SEparation) is a powder MAS heteronuclear correlation
- experiment. In the common 1H-13C implementation, molecular dynamics
- information is contained in 1H line shapes that are separated in the
- second dimension by 13C chemical shifts. Further information in:
- fid=wise(spin_system,parameters,H,R,K)
- parameters.spins -working spins, e.g. {'1H',13C'}
- parameters.hi_pwr -amplitude of high power pulses
- on the high-gamma channel, Hz
- parameters.cp_pwr -amplitude of CP pulse on each
- channel during the CP contact
- time, Hz
- parameters.cp_dur -CP contact time duration, s
