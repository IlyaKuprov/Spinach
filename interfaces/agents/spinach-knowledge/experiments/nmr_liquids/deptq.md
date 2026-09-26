# experiments/nmr_liquids/deptq.m

- Signature: `fid=deptq(spin_system,parameters,H,R,K)`

## Purpose

DEPTQ pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Syntax

```matlab
fid=deptq(spin_system,parameters,H,R,K)
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
- Note: this implementation is the DEPTQ135-style fixed first
- proton pulse variant. The beta parameter controls the
- last proton editing pulse.
- Note: use dilute.m to generate carbon isotopomers.
- Note: the sequence differs from dept.m in that quaternary carbons
- do appear.

## Implementation structure

- DEPTQ pulse sequence from:
- fid=deptq(spin_system,parameters,H,R,K)
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
