# experiments/nmr_liquids/crazed.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/crazed.m`
- Signature: `fid=crazed(spin_system,parameters,H,R,K)`
- Total lines: 138

## Purpose

CRAZED pulse sequence. Ideal analytical coherence-pathway version of the sequence described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`, `any()`.
