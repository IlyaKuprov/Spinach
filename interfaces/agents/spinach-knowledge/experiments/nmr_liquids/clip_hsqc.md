# experiments/nmr_liquids/clip_hsqc.m

- Signature: `fid=clip_hsqc(spin_system,parameters,H,R,K)`

## Purpose

CLIP-HSQC pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Syntax

```matlab
fid=clip_hsqc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. {'13C','1H'})
- parameters.J active scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay. States quadrature detection is
- used; the two components of the States signal are re-
- turned in fid.pos and fid.neg

## Implementation structure

- CLIP-HSQC pulse sequence from:
- fid=clip_hsqc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. {'13C','1H'})
- parameters.J active scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid -free induction decay. States quadrature detection is
- used; the two components of the States signal are re-
- turned in fid.pos and fid.neg
