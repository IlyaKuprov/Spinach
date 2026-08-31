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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 56-57: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 59-60: Build electron operators; implemented by `Ex=operator(spin_system,'Lx','E')`.
- Lines 64-65: Assemble the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 67-68: Loop over offsets; implemented by `dnp=zeros(size(parameters.el_offs),'like',1i)`.
- Lines 71-72: Add rotating frame offset terms to the Liouvillian; implemented by `L_curr=L+2*pi*(parameters.el_offs(n)+parameters.addshift)*Ez`.
- Lines 76-77: Add microwave irradiation along X during the 90-degree flip pulse; implemented by `L_pulse=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(0)+Ey*sind(0))`.
- Lines 79-80: Add microwave irradiation along -Y during the contact pulse; implemented by `L_contact=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(270)+Ey*sind(270))`.
- Lines 82-84: Compute and clean up NOVEL propagator; implemented by `P=propagator(spin_system,L_contact,parameters.contact_dur)* propagator(spin_system,L_pulse,parameters.pulse_dur)`.
- Lines 89-90: Add microwave irradiation along -X during the 90-degree flipback pulse; implemented by `L_flipback=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(180)+Ey*sind(180))`.
- Lines 92-93: Compute and clean up NOVEL w flipback propagator; implemented by `P=propagator(spin_system,L_flipback,parameters.pulse_dur)*P`.
- Lines 100-101: Add microwave irradiation during the contact pulse along -Y; implemented by `L_contact=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(270)+Ey*sind(270))`.
- Lines 103-104: Compute and clean up SE propagator; implemented by `P=propagator(spin_system,L_contact,parameters.contact_dur)`.
- Lines 109-110: Shot spacing delay; implemented by `P=propagator(spin_system,L_curr,parameters.shot_spacing)*P`.
- Lines 113-114: Compute the steady state; implemented by `rho=steady(spin_system,P,[],'newton')`.
- Lines 116-117: Get the observable at the steady state; implemented by `dnp(n)=gather(parameters.coil'*rho)`.

### Control flow inferred from the code

- Line 69: `for` loop over `n=1:numel(parameters.el_offs)`.
- Line 74: conditional branch on `parameters.flippulse`.
- Line 87: conditional branch on `parameters.flipback`.

### Key state/data transformations

- Lines 60: computes `Ex` using `Ex=operator(spin_system,'Lx','E')`.
- Lines 61: computes `Ey` using `Ey=operator(spin_system,'Ly','E')`.
- Lines 62: computes `Ez` using `Ez=operator(spin_system,'Lz','E')`.
- Lines 65: computes `L` using `L=H+1i*R+1i*K`.
- Lines 68: computes `dnp` using `dnp=zeros(size(parameters.el_offs),'like',1i)`.
- Lines 72: computes `L_curr` using `L_curr=L+2*pi*(parameters.el_offs(n)+parameters.addshift)*Ez`.
- Lines 77: computes `L_pulse` using `L_pulse=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(0)+Ey*sind(0))`.
- Lines 80: computes `L_contact` using `L_contact=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(270)+Ey*sind(270))`.
- Lines 83-84: computes `P` using `P=propagator(spin_system,L_contact,parameters.contact_dur)* propagator(spin_system,L_pulse,parameters.pulse_dur)`.
- Lines 90: computes `L_flipback` using `L_flipback=L_curr+2*pi*parameters.irr_powers*(Ex*cosd(180)+Ey*sind(180))`.
- Lines 114: computes `rho` using `rho=steady(spin_system,P,[],'newton')`.
- Lines 117: computes `dnp(n)` using `dnp(n)=gather(parameters.coil'*rho)`.

### Local helper functions

- Line 124: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

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
