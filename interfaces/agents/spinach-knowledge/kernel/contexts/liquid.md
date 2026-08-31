# kernel/contexts/liquid.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/liquid.m`
- Signature: `answer=liquid(spin_system,pulse_sequence,parameters,assumptions)`
- Total lines: 240

## Purpose

Liquid-phase interface to pulse sequences. Generates a Liouvillian superoperator and passes it on to the pulse sequence function, which should be supplied as a handle. This interface handles RDC mode --if the 'rdc' need is specified, it would use the order matrix supplied by the user to compute the residual anisotropies of all interactions. Syntax: answer=liquid(spin_system,pulse_sequence,parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 71-72: Show the banner; implemented by `banner(spin_system,'sequence_banner')`.
- Lines 74-75: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 77-78: Check consistency; implemented by `grumble(spin_system,pulse_sequence,parameters,assumptions)`.
- Lines 80-81: Set assumptions; implemented by `spin_system=assume(spin_system,assumptions)`.
- Lines 83-84: Get the Liouvillian; implemented by `if ismember('rdc',parameters.needs)`.
- Lines 86-87: With RDCs, first get the relaxation and kinetics; implemented by `R=relaxation(spin_system); K=kinetics(spin_system)`.
- Lines 89-90: Then do the liquid crystal averaging; implemented by `spin_system=residual(spin_system)`.
- Lines 92-93: Then get the coherent parts of the Liouvillian; implemented by `[I,Q]=hamiltonian(spin_system); H=I+orientation(Q,[0 0 0])`.
- Lines 97-98: Get the isotropic Liouvillian; implemented by `H=hamiltonian(spin_system)`.
- Lines 104-105: Get the lab frame Zeeman operator if needed; implemented by `if ismember('zeeman_op',parameters.needs)`.
- Lines 110-111: Get the thermal equilibrium if needed; implemented by `if ismember('rho_eq',parameters.needs)`.
- Lines 116-117: Process channel offsets; implemented by `H=frqoffset(spin_system,H,parameters)`.
- Lines 119-120: Get carrier operators; implemented by `C=cell(size(parameters.rframes))`.
- Lines 125-126: Apply rotating frames; implemented by `for k=1:numel(parameters.rframes)`.
- Lines 130-131: Get problem dimensions; implemented by `parameters.spc_dim=1; parameters.spn_dim=size(H,1)`.
- Lines 133-134: Report to the user; implemented by `report(spin_system,'running the pulse sequence ')`.
- Lines 136-137: Call the pulse sequence; implemented by `answer=pulse_sequence(spin_system,parameters,H,R,K)`.

### Control flow inferred from the code

- Line 84: conditional branch on `ismember('rdc',parameters.needs)`.
- Line 105: conditional branch on `ismember('zeeman_op',parameters.needs)`.
- Line 111: conditional branch on `ismember('rho_eq',parameters.needs)`.
- Line 121: `for` loop over `n=1:numel(parameters.rframes)`.
- Line 126: `for` loop over `k=1:numel(parameters.rframes)`.

### Key state/data transformations

- Lines 75: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 81: computes `spin_system` using `spin_system=assume(spin_system,assumptions)`.
- Lines 87: computes `R` using `R=relaxation(spin_system); K=kinetics(spin_system)`.
- Lines 93: computes `[I,Q]` using `[I,Q]=hamiltonian(spin_system); H=I+orientation(Q,[0 0 0])`.
- Lines 98: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 100: computes `K` using `K=kinetics(spin_system)`.
- Lines 107: computes `parameters.hzeeman` using `parameters.hzeeman=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 113: computes `parameters.rho0` using `parameters.rho0=equilibrium(spin_system)`.
- Lines 120: computes `C` using `C=cell(size(parameters.rframes))`.
- Lines 122: computes `C{n}` using `C{n}=carrier(spin_system,parameters.rframes{n}{1})`.
- Lines 131: computes `parameters.spc_dim` using `parameters.spc_dim=1; parameters.spn_dim=size(H,1)`.
- Lines 137: computes `answer` using `answer=pulse_sequence(spin_system,parameters,H,R,K)`.

### Local helper functions

- Line 142: `defaults()` — `function parameters=defaults(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'decouple')`.
  - Representative operation: `report(spin_system,'parameters.decouple field not set, assuming no decoupling.')`.
- Line 162: `grumble()` — `function grumble(spin_system,pulse_sequence,parameters,assumptions)`. Pulse sequence
  - Representative operation: `if ~isa(pulse_sequence,'function_handle')`.
  - Representative operation: `error('pulse_sequence argument must be a function handle.')`.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the
- spins that the pulse sequence works on, in
- the order of channels, e.g. {'1H','13C'}
- parameters.offset -a cell array giving
- transmitter offsets on each of the spins
- listed in parameters.spins array.
- parameters.needs -a cell array of strings specifying additional
- information required by the sequence:
- 'zeeman_op' -Zeeman part of the Hamiltonian
- in the laboratory frame, to be placed into
- parameters.hzeeman and sent to pulse sequence
- 'rdc' -triggers the processing of residual
- anisotropic couplings due to partial order
- 'rho_eq' -thermal equilibrium state at the
- specified temperature with respect to the
- isotropic part of the Hamiltonian; this is
- placed into parameters.rho0
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N,3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.* -additional subfields may be required by
- the pulse sequence -check its docs
- assumptions -context-specific assumptions ('nmr', 'epr',
- 'labframe', etc.) -see the pulse sequence
- header for information on this setting.

## Outputs

- answer -whatever it is that the pulse sequence returns.
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.

## Implementation structure

- Liquid-phase interface to pulse sequences. Generates a Liouvillian
- superoperator and passes it on to the pulse sequence function, which
- should be supplied as a handle.
- This interface handles RDC mode --if the 'rdc' need is specified,
- it would use the order matrix supplied by the user to compute the
- residual anisotropies of all interactions. Syntax:
- answer=liquid(spin_system,pulse_sequence,parameters,assumptions)
- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the
- spins that the pulse sequence works on, in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `assume()`, `ismember()`, `relaxation()`, `kinetics()`, `residual()`, `hamiltonian()`, `orientation()`, `report()`, `equilibrium()`, `frqoffset()`, `carrier()`, `rotframe()`, `pulse_sequence()`.
