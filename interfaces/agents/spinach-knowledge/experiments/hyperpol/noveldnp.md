# experiments/hyperpol/noveldnp.m

- Signature: `contact_curve=noveldnp(spin_system,parameters,H,R,K)`

## Purpose

Nuclear spin Orientation via Electron spin Locking (NOVEL) and pulsed solid effect (SE). For futher information see: Syntax (call from powder context): contact_curve=noveldnp(spin_system,parameters,H,R,K)

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
- parameters.timestep -time step of the contact curve, s
- parameters.nsteps -number of time steps in the con-
- tact curve
- parameters.flippulse -0: Solid Effect (no flip pulse)
- 1: NOVEL (90-degree flip pulse)
- Output:
- contact_curve -time dependence of the coil state

## Implementation structure

- Nuclear spin Orientation via Electron spin Locking (NOVEL) and pulsed
- solid effect (SE). For futher information see:
- Syntax (call from powder context):
- contact_curve=noveldnp(spin_system,parameters,H,R,K)
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.timestep -time step of the contact curve, s
