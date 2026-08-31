# experiments/hyperpol/xixdnp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/xixdnp.m`
- Signature: `contact_curve=xixdnp(spin_system,parameters,H,R,K)`
- Total lines: 146

## Purpose

TPPM DNP and its special case X-inverse-X (XiX) DNP experiment from (https://doi.org/10.1021/jacs.1c09900). Syntax (call from powder context): contact_curve=xixdnp(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 43-44: Build electron control operators; implemented by `Ep=operator(spin_system,'L+','E')`.
- Lines 47-48: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 50-51: Add microwave irradiation along +X; implemented by `L1=L+2*pi*parameters.irr_powers*Ex`.
- Lines 53-55: Add microwave irradiation with user-specified phase; implemented by `L2=L+2*pi*parameters.irr_powers*(Ex*cos(parameters.phase)+ Ey*sin(parameters.phase))`.
- Lines 57-58: Preallocate contact curve; implemented by `contact_curve=zeros(1,parameters.nloops+1)`.
- Lines 61-62: Grab initial condition and detection state; implemented by `rho=parameters.rho0; coil=parameters.coil`.
- Lines 64-65: Adapt to the formalism; implemented by `switch spin_system.bas.formalism`.
- Lines 67-68: Hilbert space; implemented by `case 'zeeman-hilb'`.
- Lines 70-72: Precompute complete loop propagator; implemented by `P=propagator(spin_system,L2,parameters.pulse_dur)* propagator(spin_system,L1,parameters.pulse_dur)`.
- Lines 74-75: Run and record the observable; implemented by `for k=1:parameters.nloops`.
- Lines 79-80: Liouville space; implemented by `case {'zeeman-liouv','sphten-liouv'}`.
- Lines 82-83: Precompute individual event propagators; implemented by `P1=propagator(spin_system,L1,parameters.pulse_dur)`.
- Lines 93-94: Complain and bomb out; implemented by `error('unknown formalism')`.

### Control flow inferred from the code

- Line 65: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`, `{'zeeman-liouv','sphten-liouv'}`.
- Line 75: `for` loop over `k=1:parameters.nloops`.
- Line 87: `for` loop over `k=1:parameters.nloops`.

### Key state/data transformations

- Lines 44: computes `Ep` using `Ep=operator(spin_system,'L+','E')`.
- Lines 45: computes `Ex` using `Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i`.
- Lines 48: computes `L` using `L=H+1i*R+1i*K`.
- Lines 51: computes `L1` using `L1=L+2*pi*parameters.irr_powers*Ex`.
- Lines 54-55: computes `L2` using `L2=L+2*pi*parameters.irr_powers*(Ex*cos(parameters.phase)+ Ey*sin(parameters.phase))`.
- Lines 58: computes `contact_curve` using `contact_curve=zeros(1,parameters.nloops+1)`.
- Lines 59: computes `contact_curve(1)` using `contact_curve(1)=hdot(parameters.coil,parameters.rho0)`.
- Lines 62: computes `rho` using `rho=parameters.rho0; coil=parameters.coil`.
- Lines 71-72: computes `P` using `P=propagator(spin_system,L2,parameters.pulse_dur)* propagator(spin_system,L1,parameters.pulse_dur)`.
- Lines 83: computes `P1` using `P1=propagator(spin_system,L1,parameters.pulse_dur)`.
- Lines 84: computes `P2` using `P2=propagator(spin_system,L2,parameters.pulse_dur)`.

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
- parameters.nloops -number of XiX DNP blocks
- parameters.phase -phase of the second pulse
- Output:
- contact_curve -time dependence of the coil state

## Implementation structure

- TPPM DNP and its special case X-inverse-X (XiX) DNP experiment
- from (https://doi.org/10.1021/jacs.1c09900). Syntax (call from
- powder context):
- contact_curve=xixdnp(spin_system,parameters,H,R,K)
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_dur -pulse duration, seconds

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `contact_curve()`, `hdot()`, `propagator()`, `ismatrix()`, `isfield()`, `isscalar()`.
