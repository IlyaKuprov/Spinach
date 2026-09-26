# experiments/hyperpol/beamdnp.m

- Signature: `contact_curve=beamdnp(spin_system,parameters,H,R,K)`

## Purpose

Beam DNP experiment from: Syntax (call from powder context): contact_curve=beamdnp(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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
