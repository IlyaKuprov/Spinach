# experiments/nmr_liquids/hmqc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/hmqc.m`
- Signature: `fid=hmqc(spin_system,parameters,H,R,K)`
- Total lines: 184

## Purpose

Magnitude-mode HMQC pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
fid=hmqc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths in each
- frequency direction, Hz
- parameters.npoints [F1 F2] numbers of points in
- each time direction
- parameters.spins {F1 F2} nuclei, e.g. {'15N','1H'}
- parameters.decouple_f2 nuclei to decouple in F2,
- e.g. {'15N','13C'}
- parameters.decouple_f1 nuclei that receive midpoint 180-degree
- refocusing pulses in F1, e.g. {'1H'}
- parameters.J primary scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay for magnitude mode processing
- Note: natural abundance experiments should make use of the iso-
- tope dilution functionality. See dilute.m function.

## Implementation structure

- Magnitude-mode HMQC pulse sequence from:
- fid=hmqc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths in each
- frequency direction, Hz
- parameters.npoints [F1 F2] numbers of points in
- each time direction
- parameters.spins {F1 F2} nuclei, e.g. {'15N','1H'}
- parameters.decouple_f2 nuclei to decouple in F2,
- e.g. {'15N','13C'}
- parameters.decouple_f1 nuclei that receive midpoint 180-degree
- refocusing pulses in F1, e.g. {'1H'}
- parameters.J primary scalar coupling, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `timestep()`, `decouple()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
