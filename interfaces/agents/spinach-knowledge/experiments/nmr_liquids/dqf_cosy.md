# experiments/nmr_liquids/dqf_cosy.m

- Signature: `fid=dqf_cosy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive double-quantum filtered COSY pulse sequence. Implemented as described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=dqf_cosy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos, fid.sin -components of the free induction
- decay for hypercomplex processing
- Note: the double-quantum filter is implemented as an exact
- analytical coherence-order projection, not as an explicit
- phase cycle or finite-gradient selection block.

## Implementation structure

- Phase-sensitive double-quantum filtered COSY pulse sequence.
- Implemented as described in:
- fid=dqf_cosy(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.cos, fid.sin - components of the free induction
- decay for hypercomplex processing
