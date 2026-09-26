# experiments/nmr_liquids/noesy.m

- Signature: `fid=noesy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive homonuclear NOESY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=noesy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep -sweep widths, Hz
- parameters.npoints -number of points for both dimensions
- parameters.spins -nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix -mixing time, seconds
- parameters.decouple -spins to be decoupled, specified either
- by name, e.g. {'13C','1H'}, or by a list
- of numbers, e.g. [1 2]
- parameters.rho0 -initial state; skip this and specify
- parameters.needs={'rho_eq'} to start
- from exact thermal equilibrium
- parameters.oldschool -set to 1 to disable homospoil gradient
- before the mixing time
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos,fid.sin -two components of the FID for F1 hyper-
- complex processing
- Note: this function is used for extreme simulations (proteins
- and nucleic acids) -its layout is optimised for minimum
- memory footprint rather than CPU time.
- Note: non-empty analytical decoupling is meaningful only in
- sphten-liouv formalism.

## Implementation structure

- Phase-sensitive homonuclear NOESY pulse sequence from:
- fid=noesy(spin_system,parameters,H,R,K)
- parameters.sweep -sweep widths, Hz
- parameters.npoints -number of points for both dimensions
- parameters.spins -nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix -mixing time, seconds
- parameters.decouple -spins to be decoupled, specified either
- by name, e.g. {'13C','1H'}, or by a list
- of numbers, e.g. [1 2]
- parameters.rho0 -initial state; skip this and specify
- parameters.needs={'rho_eq'} to start
