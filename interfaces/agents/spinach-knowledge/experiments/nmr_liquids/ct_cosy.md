# experiments/nmr_liquids/ct_cosy.m

- Signature: `fid=ct_cosy(spin_system,parameters,H,R,K)`

## Purpose

Constant-time COSY pulse sequence with analytical coherence selec- tion, as described in: The F1 trace is flipped at the end to preserve the conventional indirect-dimension sign for the reversed-delay implementation.

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Syntax

```matlab
fid=ct_cosy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep a two-element vector giving the
- sweep widths in F1 and F2
- parameters.npoints a two-element vector giving the
- number of points in F1 and F2
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle final pulse angle in radians, defaults
- to pi/2
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay

## Implementation structure

- Constant-time COSY pulse sequence with analytical coherence selec-
- tion, as described in:
- The F1 trace is flipped at the end to preserve the conventional
- indirect-dimension sign for the reversed-delay implementation.
- fid=ct_cosy(spin_system,parameters,H,R,K)
- parameters.sweep a two-element vector giving the
- sweep widths in F1 and F2
- parameters.npoints a two-element vector giving the
- number of points in F1 and F2
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle final pulse angle in radians, defaults
