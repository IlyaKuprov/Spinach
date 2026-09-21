# experiments/nmr_solids/fslghetcor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_solids/fslghetcor.m`
- Signature: `fid=fslghetcor(spin_system,parameters,H,R,K)`
- Total lines: 235

## Purpose

Heteronuclear correlation MAS NMR experiment with frequency-switched Lee-Goldburg homonuclear decoupling. Further details in:

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
fid=fslghetcor(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -working spins, e.g. {'1H',13C'}
- parameters.hi_pwr -amplitude of high power pulses
- on the high-gamma channel, Hz
- parameters.cp_pwr -amplitude of CP pulse on each
- channel during the CP contact
- time, Hz
- parameters.cp_dur -CP contact time duration, s
- parameters.offset -transmitter offsets on the
- two channels, Hz
- parameters.nblocks -number of FSLG blocks per
- indirect-dimension point
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.sweep -sweep width, Hz for F1, F2
- parameters.npoints -number of points in F1, F2
- H -Hamiltonian superoperator, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.sin, fid.cos -sine and cosine components
- of the States quadrature

## Implementation structure

- Heteronuclear correlation MAS NMR experiment with frequency-switched
- Lee-Goldburg homonuclear decoupling. Further details in:
- fid=fslghetcor(spin_system,parameters,H,R,K)
- parameters.spins -working spins, e.g. {'1H',13C'}
- parameters.hi_pwr -amplitude of high power pulses
- on the high-gamma channel, Hz
- parameters.cp_pwr -amplitude of CP pulse on each
- channel during the CP contact
- time, Hz
- parameters.cp_dur -CP contact time duration, s
- parameters.offset -transmitter offsets on the
- two channels, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `acos()`, `step()`, `traj_cos()`, `traj_sin()`, `ismember()`, `gpuArray()`, `clear()`, `decouple()`, `evolution()`, `dwell_times()`, `fid_cos_sin()`, `ismatrix()`, `all()`.
