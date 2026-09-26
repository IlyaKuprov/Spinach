# experiments/nmr_liquids/dept.m

- Signature: `fid=dept(spin_system,parameters,H,R,K)`

## Purpose

DEPT pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=dept(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1] Sweep width in Hz
- parameters.npoints [F1] number of points
- parameters.spins {F1,F2} nuclei, e.g. {'13C','1H'}
- parameters.J working J-coupling in Hz
- parameters.beta the angle used in the selection
- pulse, radians
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay
- Note: DEPT135 yields spectra with CH and CH3 signals in opposite phase
- to CH2 signals; DEPT90 yields spectra with only CH signals; DEPT45
- yields spectra with positive CH, CH2, and CH3 signals; quaternary
- carbons do not appear.
- Note: use dilute.m to generate carbon isotopomers.

## Implementation structure

- DEPT pulse sequence from:
- fid=dept(spin_system,parameters,H,R,K)
- parameters.sweep [F1] Sweep width in Hz
- parameters.npoints [F1] number of points
- parameters.spins {F1,F2} nuclei, e.g. {'13C','1H'}
- parameters.J working J-coupling in Hz
- parameters.beta the angle used in the selection
- pulse, radians
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid -free induction decay
