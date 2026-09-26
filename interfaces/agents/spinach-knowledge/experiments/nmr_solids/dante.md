# experiments/nmr_solids/dante.m

- Signature: `fid=dante(spin_system,parameters,H,R,K)`

## Purpose

DANTE pulse sequence. Syntax: fid=dante(spin_system,parameters,H,R,K) This function should normally be called using singlerot.m context that would provide H, R, and K.

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Parameters / inputs

- parameters.pulse_dur -duration of each pulse, seconds
- parameters.pulse_amp -amplitude of each pulse, Hz
- parameters.pulse_num -number of pulses within rotor period
- parameters.n_periods -number of rotor periods that the
- sequence is active for
- parameters.spins -working spin, specified as a
- single-element cell array
- parameters.decouple -isotopes to decouple, specified
- as a cell array
- parameters.rate -rotor frequency in Hz
- parameters.sweep -acquisition sweep width in Hz
- parameters.npoints -number of acquisition points
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+

## Outputs

- fid -free induction decay

## Implementation structure

- DANTE pulse sequence. Syntax:
- fid=dante(spin_system,parameters,H,R,K)
- This function should normally be called using singlerot.m context
- that would provide H, R, and K.
- parameters.pulse_dur -duration of each pulse, seconds
- parameters.pulse_amp -amplitude of each pulse, Hz
- parameters.pulse_num -number of pulses within rotor period
- parameters.n_periods -number of rotor periods that the
- sequence is active for
- parameters.spins -working spin, specified as a
- single-element cell array
- parameters.decouple -isotopes to decouple, specified
