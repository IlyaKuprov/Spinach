# experiments/hyperpol/xixdnp_steady.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hyperpol/xixdnp_steady.m`
- Signature: `dnp=xixdnp_steady(spin_system,parameters,H,R,K)`
- Total lines: 173

## Purpose

TPPM DNP and its special case X-inverse-X (XiX) DNP experiment from (https://doi.org/10.1021/jacs.1c09900), steady state ver- sion. Call from powder context. Syntax: dnp=xixdnp_steady(spin_system,parameters,H,R,K)

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

- Lines 55-56: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 58-59: Build electron control operators; implemented by `Ex=operator(spin_system,'Lx','E')`.
- Lines 63-64: Assemble the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 66-67: Loop over offsets; implemented by `dnp=zeros(size(parameters.el_offs),'like',1i)`.
- Lines 70-71: Add rotating frame offset terms to the Liouvillian; implemented by `L_curr=L+2*pi*(parameters.el_offs(n)+parameters.addshift)*Ez`.
- Lines 73-74: Add microwave irradiation along +X; implemented by `L1=L_curr+2*pi*parameters.irr_powers*Ex`.
- Lines 76-78: Add microwave irradiation with user-specified phase; implemented by `L2=L_curr+2*pi*parameters.irr_powers*(Ex*cos(parameters.phase)+ Ey*sin(parameters.phase))`.
- Lines 80-82: Compute the contact time block propagator; implemented by `P=propagator(spin_system,L2,parameters.pulse_dur)* propagator(spin_system,L1,parameters.pulse_dur)`.
- Lines 85-86: Send the problem to GPU if necessary; implemented by `if ismember('gpu',spin_system.sys.enable), P=gpuArray(P); end`.
- Lines 88-89: Multiply up to the loop count; implemented by `P=ppower(spin_system,P,parameters.nloops)`.
- Lines 91-92: Get the problem from GPU if necessary; implemented by `if isa(P,'gpuArray'), P=gather(P); end`.
- Lines 94-95: Add the delay after the contact time; implemented by `P=propagator(spin_system,L_curr,parameters.shot_spacing)*P`.
- Lines 98-99: Compute the stroboscopic steady state; implemented by `rho=steady(spin_system,P,[],'newton')`.
- Lines 101-102: Get the observable at the steady state; implemented by `dnp(n)=gather(parameters.coil'*rho)`.

### Control flow inferred from the code

- Line 68: `for` loop over `n=1:numel(parameters.el_offs)`.
- Line 86: conditional branch on `ismember('gpu',spin_system.sys.enable), P=gpuArray(P); end`.
- Line 92: conditional branch on `isa(P,'gpuArray'), P=gather(P); end`.

### Key state/data transformations

- Lines 59: computes `Ex` using `Ex=operator(spin_system,'Lx','E')`.
- Lines 60: computes `Ey` using `Ey=operator(spin_system,'Ly','E')`.
- Lines 61: computes `Ez` using `Ez=operator(spin_system,'Lz','E')`.
- Lines 64: computes `L` using `L=H+1i*R+1i*K`.
- Lines 67: computes `dnp` using `dnp=zeros(size(parameters.el_offs),'like',1i)`.
- Lines 71: computes `L_curr` using `L_curr=L+2*pi*(parameters.el_offs(n)+parameters.addshift)*Ez`.
- Lines 74: computes `L1` using `L1=L_curr+2*pi*parameters.irr_powers*Ex`.
- Lines 77-78: computes `L2` using `L2=L_curr+2*pi*parameters.irr_powers*(Ex*cos(parameters.phase)+ Ey*sin(parameters.phase))`.
- Lines 81-82: computes `P` using `P=propagator(spin_system,L2,parameters.pulse_dur)* propagator(spin_system,L1,parameters.pulse_dur)`.
- Lines 99: computes `rho` using `rho=steady(spin_system,P,[],'newton')`.
- Lines 102: computes `dnp(n)` using `dnp(n)=gather(parameters.coil'*rho)`.

### Local helper functions

- Line 109: `grumble()` — `function grumble(parameters,H,R,K)`.
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
- parameters.phase -phase of each second pulse in radians
- parameters.nloops -number of XiX/TPPM DNP blocks, the
- calculation is faster when this is
- an integer power of 2
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

- TPPM DNP and its special case X-inverse-X (XiX) DNP experiment
- from (https://doi.org/10.1021/jacs.1c09900), steady state ver-
- sion. Call from powder context. Syntax:
- dnp=xixdnp_steady(spin_system,parameters,H,R,K)
- H -Hamiltonian matrix, received from
- context function
- R -relaxation superoperator, received
- from context function, must be ther-
- malised to some finite temperature
- K -kinetics superoperator, received
- from context function
- parameters.irr_powers -microwave amplitude (aka electron

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `propagator()`, `clean_up()`, `ismember()`, `gpuArray()`, `ppower()`, `gather()`, `steady()`, `dnp()`, `ismatrix()`, `isfield()`, `isscalar()`.
