# experiments/nmr_solids/dante.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_solids/dante.m`
- Signature: `fid=dante(spin_system,parameters,H,R,K)`
- Total lines: 194

## Purpose

DANTE pulse sequence. Syntax: fid=dante(spin_system,parameters,H,R,K) This function should normally be called using singlerot.m context that would provide H, R, and K.

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 48-49: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 51-52: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 54-55: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 58-59: Apply the decoupling; implemented by `[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 61-62: Timing parameters; implemented by `rotor_period=abs(1/parameters.rate)`.
- Lines 65-66: Bomb out if the schedule makes no sense; implemented by `if (cycle_length-parameters.pulse_dur)<0`.
- Lines 70-71: Get the state vector going; implemented by `rho=parameters.rho0`.
- Lines 73-74: Loop over rotor cycles; implemented by `for k=1:parameters.n_periods`.
- Lines 76-77: Loop over pulses within the rotor cycle; implemented by `for n=1:parameters.pulse_num`.
- Lines 79-81: Pulse; implemented by `rho=evolution(spin_system,L+2*pi*parameters.pulse_amp*Lx, [],rho,parameters.pulse_dur,1,'final')`.
- Lines 83-85: Free evolution; implemented by `rho=evolution(spin_system,L, [],rho,cycle_length-parameters.pulse_dur,1,'final')`.
- Lines 91-93: Run the evolution and watch the coil state; implemented by `fid=evolution(spin_system,L,parameters.coil,rho, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Control flow inferred from the code

- Line 66: conditional branch on `(cycle_length-parameters.pulse_dur)<0`.
- Line 74: `for` loop over `k=1:parameters.n_periods`.
- Line 77: `for` loop over `n=1:parameters.pulse_num`.

### Key state/data transformations

- Lines 52: computes `L` using `L=H+1i*R+1i*K`.
- Lines 55: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 56: computes `Lx` using `Lx=(Lp+Lp')/2; Lx=kron(speye(parameters.spc_dim),Lx)`.
- Lines 59: computes `[L,parameters.rho0]` using `[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 62: computes `rotor_period` using `rotor_period=abs(1/parameters.rate)`.
- Lines 63: computes `cycle_length` using `cycle_length=rotor_period/parameters.pulse_num`.
- Lines 71: computes `rho` using `rho=parameters.rho0`.
- Lines 92-93: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 98: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.pulse_dur -duration of each pulse, seconds
- parameters.pulse_amp -amplitude of each pulse, Hz
- parameters.pulse_num -number of pulses within rotor period
- parameters.n_periods -number of rotor periods that the
- sequence is active for
- parameters.spins -working spin, specified as a
- single-element cell array
- parameters.decouple -isotopes to decouple, specified
- as a cell array
- parameters.rate -rotor frequency in Hz
- parameters.sweep -acquisition sweep width in Hz
- parameters.npoints -number of acquisition points
- parameters.spc_dim -Fokker-Planck spatial dimension
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+

## Outputs

- fid -free induction decay

## Implementation structure

- DANTE pulse sequence. Syntax:
- fid=dante(spin_system,parameters,H,R,K)
- This function should normally be called using singlerot.m context
- that would provide H, R, and K.
- parameters.pulse_dur -duration of each pulse, seconds
- parameters.pulse_amp -amplitude of each pulse, Hz
- parameters.pulse_num -number of pulses within rotor period
- parameters.n_periods -number of rotor periods that the
- sequence is active for
- parameters.spins -working spin, specified as a
- single-element cell array
- parameters.decouple -isotopes to decouple, specified

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `decouple()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ismember()`, `spins()`, `any()`, `cellfun()`.
