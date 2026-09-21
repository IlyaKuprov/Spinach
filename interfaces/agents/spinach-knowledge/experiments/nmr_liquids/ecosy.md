# experiments/nmr_liquids/ecosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/ecosy.m`
- Signature: `fid=ecosy(spin_system,parameters,H,R,K)`
- Total lines: 128

## Purpose

Phase-sensitive E.COSY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
fid=ecosy(spin_system,parameters,H,R,K)
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

- fid.cos, fid.sin -real and imaginary components of the
- States quadrature signal
- Notes: implemented up to six-quantum orders as per the original
- paper. Let us know if you need more.

## Implementation structure

- Phase-sensitive E.COSY pulse sequence from:
- fid=ecosy(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.cos, fid.sin -real and imaginary components of the
- States quadrature signal
- paper. Let us know if you need more.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `evolution()`, `step()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`, `any()`.
