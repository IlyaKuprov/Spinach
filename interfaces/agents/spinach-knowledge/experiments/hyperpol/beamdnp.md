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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 43-44: Build electron control operators; implemented by `Ep=operator(spin_system,'L+','E')`.
- Lines 47-48: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 50-51: Microwave irradiation along Y during the 90-degree flip pulse; implemented by `L_pulse=L+2*pi*parameters.irr_powers*Ey`.
- Lines 53-54: Calculate the length of the 90-degree flip pulse; implemented by `pulse_dur=1/(4*parameters.irr_powers)`.
- Lines 56-57: Apply the 90-degree flip pulse; implemented by `rho=evolution(spin_system,L_pulse,[],parameters.rho0,pulse_dur,1,'final')`.
- Lines 59-60: Make repeated pulse Liouvillians; implemented by `L1=L+2*pi*parameters.irr_powers*Ex`.
- Lines 63-64: Preallocate contact curve; implemented by `contact_curve=zeros(1,parameters.nloops+1)`.
- Lines 67-68: Adapt to the formalism; implemented by `switch spin_system.bas.formalism`.
- Lines 70-71: Hilbert space; implemented by `case 'zeeman-hilb'`.
- Lines 73-75: Precompute the complete loop propagator; implemented by `P=propagator(spin_system,L2,parameters.pulse_dur(2))* propagator(spin_system,L1,parameters.pulse_dur(1))`.
- Lines 77-78: Run the sequence; implemented by `for k=1:parameters.nloops`.
- Lines 80-81: Run a loop and record the observable; implemented by `rho=P*rho*P'; contact_curve(1+k)=hdot(parameters.coil,rho)`.
- Lines 85-86: Liouville space; implemented by `case {'zeeman-liouv','sphten-liouv'}`.
- Lines 88-89: Precompute individual event propagators; implemented by `P1=propagator(spin_system,L1,parameters.pulse_dur(1))`.
- Lines 95-96: Run a loop and record the observable; implemented by `rho=P2*(P1*rho); contact_curve(1+k)=hdot(parameters.coil,rho)`.
- Lines 102-103: Complain and bomb out; implemented by `error('unknown formalism.')`.

### Control flow inferred from the code

- Line 68: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`, `{'zeeman-liouv','sphten-liouv'}`.
- Line 78: `for` loop over `k=1:parameters.nloops`.
- Line 93: `for` loop over `k=1:parameters.nloops`.

### Key state/data transformations

- Lines 44: computes `Ep` using `Ep=operator(spin_system,'L+','E')`.
- Lines 45: computes `Ex` using `Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i`.
- Lines 48: computes `L` using `L=H+1i*R+1i*K`.
- Lines 51: computes `L_pulse` using `L_pulse=L+2*pi*parameters.irr_powers*Ey`.
- Lines 54: computes `pulse_dur` using `pulse_dur=1/(4*parameters.irr_powers)`.
- Lines 57: computes `rho` using `rho=evolution(spin_system,L_pulse,[],parameters.rho0,pulse_dur,1,'final')`.
- Lines 60: computes `L1` using `L1=L+2*pi*parameters.irr_powers*Ex`.
- Lines 61: computes `L2` using `L2=L-2*pi*parameters.irr_powers*Ex`.
- Lines 64: computes `contact_curve` using `contact_curve=zeros(1,parameters.nloops+1)`.
- Lines 65: computes `contact_curve(1)` using `contact_curve(1)=hdot(parameters.coil,parameters.rho0)`.
- Lines 74-75: computes `P` using `P=propagator(spin_system,L2,parameters.pulse_dur(2))* propagator(spin_system,L1,parameters.pulse_dur(1))`.
- Lines 89: computes `P1` using `P1=propagator(spin_system,L1,parameters.pulse_dur(1))`.
- Lines 90: computes `P2` using `P2=propagator(spin_system,L2,parameters.pulse_dur(2))`.

### Local helper functions

- Line 110: `grumble()` — `function grumble(parameters,H,R,K)`.
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
