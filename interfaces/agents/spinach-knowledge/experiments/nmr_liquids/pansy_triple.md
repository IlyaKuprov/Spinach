# experiments/nmr_liquids/pansy_triple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/pansy_triple.m`
- Signature: `fid=pansy_triple(spin_system,parameters,H,R,K)`
- Total lines: 149

## Purpose

Triple-channel PANSY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
fid=pansy_triple(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -three nuclei on which the sequence
- runs, e.g. {'1H','13C','15N'}
- parameters.sweep -a vector with three sweep widths
- in Hz
- parameters.npoints -a vector of three integers specify-
- ing point count in each dimension
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.aa -COSY FID on F1 nuclei
- fid.ab -COSY FID on F1,F2 nuclei
- fid.ac -COSY FID on F1,F3 nuclei
- Note: decoupling with respect to any of the working nuclei is
- impossible in this pulse sequence.
- Note: this is the magnitude-mode analytical-pathway version.
- Gradient echo/anti-echo PANSY variants should be implemented
- as separate pulse sequence functions.
- Andrew Porter

## Implementation structure

- Triple-channel PANSY pulse sequence from:
- fid=pansy_triple(spin_system,parameters,H,R,K)
- parameters.spins -three nuclei on which the sequence
- runs, e.g. {'1H','13C','15N'}
- parameters.sweep -a vector with three sweep widths
- in Hz
- parameters.npoints -a vector of three integers specify-
- ing point count in each dimension
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.aa -COSY FID on F1 nuclei

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `coherence()`, `evolution()`, `timesteps()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
