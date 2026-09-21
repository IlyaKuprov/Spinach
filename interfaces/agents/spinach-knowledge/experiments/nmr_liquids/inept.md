# experiments/nmr_liquids/inept.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/inept.m`
- Signature: `fid=inept(spin_system,parameters,H,R,K)`
- Total lines: 146

## Purpose

Non-refocused INEPT pulse sequence. This returns the directly acquired coupled antiphase spectrum rather than a refocused, broadband-decoupled INEPT variant. Implemented as here:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
fid=inept(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1] Sweep width in Hz
- parameters.npoints [F1] number of points
- parameters.spins {F1 F2} working nuclei,
- e.g. {'15N','1H'}
- parameters.J working scalar coupling
- in Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay
- Note: use dilute() to generate carbon isotopomers.
- Andrew Porter, Ilya Kuprov

## Implementation structure

- Non-refocused INEPT pulse sequence. This returns the directly
- acquired coupled antiphase spectrum rather than a refocused,
- broadband-decoupled INEPT variant. Implemented as here:
- fid=inept(spin_system,parameters,H,R,K)
- parameters.sweep [F1] Sweep width in Hz
- parameters.npoints [F1] number of points
- parameters.spins {F1 F2} working nuclei,
- e.g. {'15N','1H'}
- parameters.J working scalar coupling
- in Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `equilibrium()`, `state()`, `operator()`, `step()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`, `strcmp()`, `any()`.
