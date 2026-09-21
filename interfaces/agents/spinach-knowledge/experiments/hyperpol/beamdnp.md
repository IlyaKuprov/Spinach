# experiments/hyperpol/beamdnp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/beamdnp.m`
- Signature: `contact_curve=beamdnp(spin_system,parameters,H,R,K)`
- Total lines: 151

## Purpose

Beam DNP experiment from: Syntax (call from powder context): contact_curve=beamdnp(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_dur -two pulse durations, seconds
- parameters.nloops -number of BEAM DNP blocks
- Output:
- contact_curve -time dependence of the coil state

## Implementation structure

- Beam DNP experiment from:
- Syntax (call from powder context):
- contact_curve=beamdnp(spin_system,parameters,H,R,K)
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_dur -two pulse durations, seconds
- parameters.nloops -number of BEAM DNP blocks

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `evolution()`, `contact_curve()`, `hdot()`, `propagator()`, `ismatrix()`, `isfield()`, `isscalar()`, `isrow()`, `any()`.
