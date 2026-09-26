# experiments/nmr_liquids/crazed.m

- Signature: `fid=crazed(spin_system,parameters,H,R,K)`

## Purpose

CRAZED pulse sequence. Ideal analytical coherence-pathway version of the sequence described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=crazed(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle
- parameters.rho0 initial condition
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay
- Note: the gradient selection is represented by explicit coherence
- projection onto the +2 double-quantum branch followed by the
- +1 observable single-quantum branch.

## Implementation structure

- CRAZED pulse sequence. Ideal analytical coherence-pathway version
- of the sequence described in:
- fid=crazed(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle
- parameters.rho0 initial condition
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
