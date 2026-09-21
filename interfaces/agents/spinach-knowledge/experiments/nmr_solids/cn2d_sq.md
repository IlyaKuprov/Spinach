# experiments/nmr_solids/cn2d_sq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_solids/cn2d_sq.m`
- Signature: `fid=cn2d_sq(spin_system,parameters,H,R,K)`
- Total lines: 141

## Purpose

Single-quantum version of the 13C-detected 14N-13C MAS 2D correlation experiment described by Jarvis, Haies, Williamson and Carravetta in

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
fid=cn2d_sq(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -isotopes to which the sequence is
- applied, specified as a cell array
- with 14N first, and 13C second
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.sweep -sweep widths in the two dimensions, Hz
- parameters.npoints -numbers of points in the two dimensions
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.rf_pwr -RF power on 14N, Hz
- parameters.rf_dur -RF pulse duration on 14N, seconds

## Outputs

- fid.sin
- fid.cos -sine and cosine components
- of the States quadrature

## Implementation structure

- Single-quantum version of the 13C-detected 14N-13C MAS 2D correlation
- experiment described by Jarvis, Haies, Williamson and Carravetta in
- fid=cn2d_sq(spin_system,parameters,H,R,K)
- parameters.spins -isotopes to which the sequence is
- applied, specified as a cell array
- with 14N first, and 13C second
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.sweep -sweep widths in the two dimensions, Hz
- parameters.npoints -numbers of points in the two dimensions
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.rf_pwr -RF power on 14N, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `evolution()`, `coherence()`, `timestep()`, `step()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`.
