# experiments/hyperpol/topdnp_steady.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/topdnp_steady.m`
- Signature: `dnp=topdnp_steady(spin_system,parameters,H,R,K)`
- Total lines: 170

## Purpose

Time-optimised pulsed DNP experiment from: (a steady state version). Syntax (call from powder context): dnp=topdnp_steady(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Hyperpolarisation experiment implementations. They propagate driven electron-nuclear systems under microwave irradiation, MAS, relaxation, and repetition until transient or steady-state observables are assembled.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 56-57: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 59-60: Build electron operators; implemented by `Ex=operator(spin_system,'Lx','E')`.
- Lines 63-64: Assemble the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 66-67: Loop over offsets; implemented by `dnp=zeros(size(parameters.el_offs),'like',1i)`.
- Lines 70-71: Add rotating frame offset terms to the Liouvillian; implemented by `L_curr=L+2*pi*(parameters.el_offs(n)+parameters.addshift)*Ez`.
- Lines 73-74: Add microwave irradiation along +X; implemented by `L1=L_curr+2*pi*parameters.irr_powers*Ex`.
- Lines 76-78: Precompute and clean up TOP DNP propagator; implemented by `P=propagator(spin_system,L_curr,parameters.delay_dur)* propagator(spin_system,L1,parameters.pulse_dur)`.
- Lines 81-82: Send the problem to GPU if necessary; implemented by `if ismember('gpu',spin_system.sys.enable), P=gpuArray(P); end`.
- Lines 84-85: Multiply up to the loop count; implemented by `P=ppower(spin_system,P,parameters.nloops)`.
- Lines 87-88: Get the problem from GPU if necessary; implemented by `if isa(P,'gpuArray'), P=gather(P); end`.
- Lines 90-91: Shot spacing delay; implemented by `P=propagator(spin_system,L_curr,parameters.shot_spacing)*P`.
- Lines 94-95: Compute the steady state; implemented by `rho=steady(spin_system,P,[],'newton')`.
- Lines 97-98: Get the observable at the steady state; implemented by `dnp(n)=gather(parameters.coil'*rho)`.

### Control flow inferred from the code

- Line 68: `for` loop over `n=1:numel(parameters.el_offs)`.
- Line 82: conditional branch on `ismember('gpu',spin_system.sys.enable), P=gpuArray(P); end`.
- Line 88: conditional branch on `isa(P,'gpuArray'), P=gather(P); end`.

### Key state/data transformations

- Lines 60: computes `Ex` using `Ex=operator(spin_system,'Lx','E')`.
- Lines 61: computes `Ez` using `Ez=operator(spin_system,'Lz','E')`.
- Lines 64: computes `L` using `L=H+1i*R+1i*K`.
- Lines 67: computes `dnp` using `dnp=zeros(size(parameters.el_offs),'like',1i)`.
- Lines 71: computes `L_curr` using `L_curr=L+2*pi*(parameters.el_offs(n)+parameters.addshift)*Ez`.
- Lines 74: computes `L1` using `L1=L_curr+2*pi*parameters.irr_powers*Ex`.
- Lines 77-78: computes `P` using `P=propagator(spin_system,L_curr,parameters.delay_dur)* propagator(spin_system,L1,parameters.pulse_dur)`.
- Lines 95: computes `rho` using `rho=steady(spin_system,P,[],'newton')`.
- Lines 98: computes `dnp(n)` using `dnp(n)=gather(parameters.coil'*rho)`.

### Local helper functions

- Line 105: `grumble()` — `function grumble(parameters,H,R,K)`.
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
- parameters.coil -detection state vector
- parameters.pulse_dur -pulse duration, seconds
- parameters.delay_dur -delay duration, seconds
- parameters.nloops -number of XiX/TPPM DNP blocks,
- must be an integer power of 2
- parameters.shot_spacing -delay between microwave irradiation
- periods, seconds
- parameters.addshift -shift of the centre of the field
- profile, Hz
- parameters.el_offs -microwave resonance offsets, a vector
- of frequencies in Hz
- Output:
- dnp -steady state observable on the de-
- tection state vector as a function
- of microwave resonance offset

## Implementation structure

- Time-optimised pulsed DNP experiment from:
- (a steady state version). Syntax (call from powder context):
- dnp=topdnp_steady(spin_system,parameters,H,R,K)
- H -Hamiltonian matrix, received from
- context function
- R -relaxation superoperator, received
- from context function, must be ther-
- malised to some finite temperature
- K -kinetics superoperator, received
- from context function
- parameters.irr_powers -microwave amplitude (aka electron
- nutation frequency), Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `propagator()`, `clean_up()`, `ismember()`, `gpuArray()`, `ppower()`, `gather()`, `steady()`, `dnp()`, `ismatrix()`, `isfield()`, `isscalar()`.
