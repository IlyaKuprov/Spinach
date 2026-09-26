# experiments/nmr_liquids/inadequate_2d.m

- Signature: `fid=inadequate_2d(spin_system,parameters,H,R,K)`

## Purpose

2D INADEQUATE pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Syntax

```matlab
fid=inadequate_2d(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins active nucleus, e.g. {'13C'}
- parameters.decouple optional nuclei to decouple, e.g. {'1H'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos,fid.sin -two components of the States signal
- Notes: use dilute.m to generate carbon pair isotopomers. The F1
- axis is a double-quantum frequency coordinate.
- Theresa Hune
- Christian Griesinger

## Implementation structure

- 2D INADEQUATE pulse sequence from:
- fid=inadequate_2d(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins active nucleus, e.g. {'13C'}
- parameters.decouple optional nuclei to decouple, e.g. {'1H'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.cos,fid.sin - two components of the States signal
- axis is a double-quantum frequency coordinate.
