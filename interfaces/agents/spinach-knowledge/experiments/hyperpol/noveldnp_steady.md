# experiments/hyperpol/noveldnp_steady.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/noveldnp_steady.m`
- Signature: `dnp=noveldnp_steady(spin_system,parameters,H,R,K)`
- Total lines: 193

## Purpose

Nuclear spin Orientation via Electron spin Locking (NOVEL) and pulsed solid effect (SE), steady-state version. For futher information see: Syntax (call from powder context): dnp=noveldnp_steady(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- H -Hamiltonian matrix, received from
- context function
- R -relaxation superoperator, received
- from context function, must be ther-
- malised to some finite temperature
- K -kinetics superoperator, received
- from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz
- parameters.coil -detection state(s)
- parameters.contact_dur -contact time, seconds
- parameters.shot_spacing -delay between microwave irradiation periods
- parameters.flippulse -0: Solid Effect (no flip pulse)
- 1: NOVEL (90-degree flip pulse)
- parameters.flipback -0: NOVEL without flipback pulse
- 1: NOVEL with flipback pulse
- parameters.addshift -shift to center the field profile
- parameters.el_offs -microwave resonance offsets
- Output:
- dnp -steady state observable on the de-
- tection state vector as a function
- of microwave resonance offset

## Implementation structure

- Nuclear spin Orientation via Electron spin Locking (NOVEL) and pulsed
- solid effect (SE), steady-state version. For futher information see:
- Syntax (call from powder context):
- dnp=noveldnp_steady(spin_system,parameters,H,R,K)
- H -Hamiltonian matrix, received from
- context function
- R -relaxation superoperator, received
- from context function, must be ther-
- malised to some finite temperature
- K -kinetics superoperator, received
- from context function
- parameters.irr_powers -microwave amplitude (aka electron

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `cosd()`, `sind()`, `propagator()`, `clean_up()`, `steady()`, `dnp()`, `gather()`, `ismatrix()`, `isfield()`, `isscalar()`.
