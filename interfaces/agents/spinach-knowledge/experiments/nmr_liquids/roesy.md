# experiments/nmr_liquids/roesy.m

- Signature: `fid=roesy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive homonuclear ROESY pulse sequence, assuming ideal spin-lock, described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=roesy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep -a vector with sweep widths
- in F1 and F2 directions, Hz
- parameters.npoints -a vector with point count
- in F1 and F2 directions
- parameters.spins -nuclei on which the sequence
- runs, e.g. {'1H'}
- parameters.tmix -mixing time, seconds
- parameters.rho0 -initial state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos, fid.sin -components of the free induction
- decay for hypercomplex processing
- Note: this ideal spin-lock model does not represent finite RF
- amplitude, RF offset, Hartmann-Hahn matching errors, or
- explicit RF phase transients.

## Implementation structure

- Phase-sensitive homonuclear ROESY pulse sequence, assuming ideal
- spin-lock, described in:
- fid=roesy(spin_system,parameters,H,R,K)
- parameters.sweep -a vector with sweep widths
- in F1 and F2 directions, Hz
- parameters.npoints -a vector with point count
- in F1 and F2 directions
- parameters.spins -nuclei on which the sequence
- runs, e.g. {'1H'}
- parameters.tmix -mixing time, seconds
- parameters.rho0 -initial state
- H -Hamiltonian matrix, received from context function
