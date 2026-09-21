# experiments/overtone/overtone_dante.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/overtone/overtone_dante.m`
- Signature: `spectrum=overtone_dante(spin_system,parameters,H,R,K)`
- Total lines: 189

## Purpose

Overtone DANTE experiment with frequency-domain acquisition.

## Physical / mathematical content

- Overtone experiment implementations. These routines excite or detect high-order quadrupolar transitions and therefore combine non-secular quadrupolar terms, MAS or field effects, and specialised detection pathways.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
spectrum=overtone_dante(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.pulse_dur -duration of the pulse, seconds
- parameters.pulse_amp -amplitude of the pulse, rad/s
- parameters.pulse_num -number of pulses within the
- rotor period
- parameters.n_periods -number of rotor periods that the
- sequence is active for
- parameters.spins -overtone-active nucleus, specified
- as a single-element cell array
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.Lx -X Zeeman operator on the
- quadrupolar nucleus
- parameters.rf_frq -pulse frequency offset from the
- overtone frequency, Hz
- parameters.rate -rotor frequency in Hz
- parameters.sweep -acquisition sweep range, Hz
- parameters.npoints -number of acquisition points
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -overtone spectrum

## Implementation structure

- Overtone DANTE experiment with frequency-domain acquisition.
- spectrum=overtone_dante(spin_system,parameters,H,R,K)
- parameters.pulse_dur -duration of the pulse, seconds
- parameters.pulse_amp -amplitude of the pulse, rad/s
- parameters.pulse_num -number of pulses within the
- rotor period
- parameters.n_periods -number of rotor periods that the
- sequence is active for
- parameters.spins -overtone-active nucleus, specified
- as a single-element cell array
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.Lx -X Zeeman operator on the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `speye()`, `average()`, `propagator()`, `clean_up()`, `multiprop()`, `overtone_a()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`.
