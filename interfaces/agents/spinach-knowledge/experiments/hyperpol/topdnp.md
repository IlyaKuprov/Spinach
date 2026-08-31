# experiments/hyperpol/topdnp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/topdnp.m`
- Signature: `contact_curve=topdnp(spin_system,parameters,H,R,K)`
- Total lines: 163

## Purpose

Time-optimised pulsed DNP experiment from: Syntax (call from powder context): contact_curve=topdnp(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 45-46: Build electron control operators; implemented by `Ep=operator(spin_system,'L+','E'); Ex=(Ep+Ep')/2`.
- Lines 48-49: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 51-52: Add microwave irradiation along X; implemented by `L_pulse=L+2*pi*parameters.irr_powers*Ex`.
- Lines 54-55: Preallocate contact curve; implemented by `contact_curve=zeros(1,parameters.nloops+1)`.
- Lines 58-59: Adapt to the formalism; implemented by `switch spin_system.bas.formalism`.
- Lines 61-62: Hilbert space; implemented by `case 'zeeman-hilb'`.
- Lines 64-66: Precompute TOP DNP loop propagator; implemented by `P=propagator(spin_system,L,parameters.delay_dur)* propagator(spin_system,L_pulse,parameters.pulse_dur)`.
- Lines 68-69: Run TOP DNP sequence; implemented by `rho=parameters.rho0`.
- Lines 72-73: Run a loop and record the observable; implemented by `rho=P*rho*P'; contact_curve(1+k)=hdot(parameters.coil,rho)`.
- Lines 77-78: Liouville space; implemented by `case {'zeeman-liouv','sphten-liouv'}`.
- Lines 84-85: Run a loop and record the observable; implemented by `rho=step(spin_system,L_pulse,rho,parameters.pulse_dur)`.
- Lines 93-94: Complain and bomb out; implemented by `error('unknown formalism')`.

### Control flow inferred from the code

- Line 59: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`, `{'zeeman-liouv','sphten-liouv'}`.
- Line 70: `for` loop over `k=1:parameters.nloops`.
- Line 82: `for` loop over `k=1:parameters.nloops`.

### Key state/data transformations

- Lines 46: computes `Ep` using `Ep=operator(spin_system,'L+','E'); Ex=(Ep+Ep')/2`.
- Lines 49: computes `L` using `L=H+1i*R+1i*K`.
- Lines 52: computes `L_pulse` using `L_pulse=L+2*pi*parameters.irr_powers*Ex`.
- Lines 55: computes `contact_curve` using `contact_curve=zeros(1,parameters.nloops+1)`.
- Lines 56: computes `contact_curve(1)` using `contact_curve(1)=hdot(parameters.coil,parameters.rho0)`.
- Lines 65-66: computes `P` using `P=propagator(spin_system,L,parameters.delay_dur)* propagator(spin_system,L_pulse,parameters.pulse_dur)`.
- Lines 69: computes `rho` using `rho=parameters.rho0`.
- Lines 87: computes `contact_curve(1+k)` using `contact_curve(1+k)=hdot(parameters.coil,rho)`.

### Local helper functions

- Line 101: `grumble()` — `function grumble(parameters,H,R,K)`.
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
- parameters.pulse_dur -pulse duration, seconds
- parameters.delay_dur -delay_duration, seconds
- parameters.nloops -number of TOP DNP loops
- Output:
- contact_curve -time dependence of the coil state

## Implementation structure

- Time-optimised pulsed DNP experiment from:
- Syntax (call from powder context):
- contact_curve=topdnp(spin_system,parameters,H,R,K)
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_dur -pulse duration, seconds
- parameters.delay_dur -delay_duration, seconds

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `contact_curve()`, `hdot()`, `propagator()`, `step()`, `ismatrix()`, `isfield()`, `isscalar()`.
