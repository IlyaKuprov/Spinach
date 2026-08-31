# experiments/hyperpol/noveldnp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/noveldnp.m`
- Signature: `contact_curve=noveldnp(spin_system,parameters,H,R,K)`
- Total lines: 147

## Purpose

Nuclear spin Orientation via Electron spin Locking (NOVEL) and pulsed solid effect (SE). For futher information see: Syntax (call from powder context): contact_curve=noveldnp(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 49-50: Build electron control operators; implemented by `Ep=operator(spin_system,'L+','E')`.
- Lines 53-54: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 56-57: 90-degree flip pulse; implemented by `if parameters.flippulse`.
- Lines 59-60: Add microwave irradiation along X during the 90-degree flip pulse; implemented by `L_pulse=L+2*pi*parameters.irr_powers*(Ex*cosd(0)+Ey*sind(0))`.
- Lines 62-64: Apply the 90-degree flip pulse; implemented by `rho=evolution(spin_system,L_pulse,[],parameters.rho0, parameters.pulse_dur,1,'final')`.
- Lines 68-69: Do nothing; implemented by `rho=parameters.rho0`.
- Lines 73-74: Set the spin lock direction to -Y; implemented by `L_slock=L+2*pi*parameters.irr_powers*(Ex*cosd(270)+Ey*sind(270))`.
- Lines 76-78: Compute the contact curve; implemented by `contact_curve=evolution(spin_system,L_slock,parameters.coil,rho, parameters.timestep,parameters.nsteps,'observable')`.

### Control flow inferred from the code

- Line 57: conditional branch on `parameters.flippulse`.

### Key state/data transformations

- Lines 50: computes `Ep` using `Ep=operator(spin_system,'L+','E')`.
- Lines 51: computes `Ex` using `Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i`.
- Lines 54: computes `L` using `L=H+1i*R+1i*K`.
- Lines 60: computes `L_pulse` using `L_pulse=L+2*pi*parameters.irr_powers*(Ex*cosd(0)+Ey*sind(0))`.
- Lines 63-64: computes `rho` using `rho=evolution(spin_system,L_pulse,[],parameters.rho0, parameters.pulse_dur,1,'final')`.
- Lines 74: computes `L_slock` using `L_slock=L+2*pi*parameters.irr_powers*(Ex*cosd(270)+Ey*sind(270))`.
- Lines 77-78: computes `contact_curve` using `contact_curve=evolution(spin_system,L_slock,parameters.coil,rho, parameters.timestep,parameters.nsteps,'observable')`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `cosd()`, `sind()`, `evolution()`, `ismatrix()`, `isfield()`, `isscalar()`.
