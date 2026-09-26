# experiments/nmr_liquids/hoesy.m

- Signature: `fid=hoesy(spin_system,parameters,H,R,K)`

## Purpose

Phase-sensitive heteronuclear NOESY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=hoesy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep two sweep widths, Hz
- parameters.npoints number of FID points for both
- dimensions
- parameters.spins nuclei on which the sequence runs,
- e.g. {'15N','13C'}
- parameters.decouple_f1 nuclei to decouple in F1, e.g.
- {'1H','13C'}
- parameters.tmix mixing time, seconds
- parameters.rho0 initial state
- parameters.needs should be set to {'rho_eq'}, this
- sequence needs the thermal equili-
- brium state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos,fid.sin -two components of the FID for F1 hyper-
- complex processing
- Note: this is an ideal heteronuclear NOESY model. Gradient and
- diffusion attenuation, finite-pulse losses, and experimental
- normalisation are outside this pulse sequence function.

## Implementation structure

- Phase-sensitive heteronuclear NOESY pulse sequence from:
- fid=hoesy(spin_system,parameters,H,R,K)
- parameters.sweep two sweep widths, Hz
- parameters.npoints number of FID points for both
- dimensions
- parameters.spins nuclei on which the sequence runs,
- e.g. {'15N','13C'}
- parameters.decouple_f1 nuclei to decouple in F1, e.g.
- {'1H','13C'}
- parameters.tmix mixing time, seconds
- parameters.rho0 initial state
- parameters.needs should be set to {'rho_eq'}, this
