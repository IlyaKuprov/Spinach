# kernel/contexts/crystal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/crystal.m`
- Signature: `answer=crystal(spin_system,pulse_sequence,parameters,assumptions)`
- Total lines: 259

## Purpose

Single-crystal interface to pulse sequences. Generates a Liouvillian superoperator and passes it on to the pulse sequence function, which should be supplied as a handle. Syntax: answer=crystal(spin_system,pulse_sequence,parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 76-77: Show the banner; implemented by `banner(spin_system,'sequence_banner')`.
- Lines 79-80: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 82-83: Check consistency; implemented by `grumble(spin_system,pulse_sequence,parameters,assumptions)`.
- Lines 85-86: Report to the user; implemented by `report(spin_system,'building the Liouvillian ')`.
- Lines 88-89: Set the assumptions; implemented by `spin_system=assume(spin_system,assumptions)`.
- Lines 91-92: Get the rotational expansion; implemented by `[I,Q]=hamiltonian(spin_system)`.
- Lines 94-95: Get the lab frame Zeeman operator if needed; implemented by `if ismember('zeeman_op',parameters.needs)`.
- Lines 100-101: Compute the thermal equilibrium at the current orientation; implemented by `if ismember('aniso_eq',parameters.needs)`.
- Lines 106-107: Apply the offsets; implemented by `I=frqoffset(spin_system,I,parameters)`.
- Lines 109-110: Build the Hamiltonian and tidy up; implemented by `H=I+orientation(Q,parameters.orientation)`.
- Lines 113-114: Get the lab frame Zeeman operator at the current orientation; implemented by `if ismember('zeeman_op',parameters.needs)`.
- Lines 119-120: Apply rotating frames; implemented by `for k=1:numel(parameters.rframes)`.
- Lines 122-123: Get the carrier operator; implemented by `C=carrier(spin_system,parameters.rframes{k}{1})`.
- Lines 125-127: Compute the rotating frame transformation; implemented by `H=rotframe(spin_system,C,H,parameters.rframes{k}{1}, parameters.rframes{k}{2})`.
- Lines 131-132: Build relaxation at the orientation specified; implemented by `R=relaxation(spin_system,parameters.orientation)`.
- Lines 134-135: Build kinetics; implemented by `K=kinetics(spin_system)`.
- Lines 137-138: Get problem dimensions; implemented by `parameters.spc_dim=1; parameters.spn_dim=size(H,1)`.
- Lines 140-141: Report to the user; implemented by `report(spin_system,'running the pulse sequence ')`.

### Control flow inferred from the code

- Line 95: conditional branch on `ismember('zeeman_op',parameters.needs)`.
- Line 101: conditional branch on `ismember('aniso_eq',parameters.needs)`.
- Line 114: conditional branch on `ismember('zeeman_op',parameters.needs)`.
- Line 120: `for` loop over `k=1:numel(parameters.rframes)`.

### Key state/data transformations

- Lines 80: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 89: computes `spin_system` using `spin_system=assume(spin_system,assumptions)`.
- Lines 92: computes `[I,Q]` using `[I,Q]=hamiltonian(spin_system)`.
- Lines 97: computes `[ZI,ZQ]` using `[ZI,ZQ]=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 102: computes `[HL,QL]` using `[HL,QL]=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 103: computes `parameters.rho0` using `parameters.rho0=equilibrium(spin_system,HL,QL,parameters.orientation)`.
- Lines 107: computes `I` using `I=frqoffset(spin_system,I,parameters)`.
- Lines 110: computes `H` using `H=I+orientation(Q,parameters.orientation)`.
- Lines 115: computes `Z` using `Z=ZI+orientation(ZQ,parameters.orientation)`.
- Lines 116: computes `parameters.hzeeman` using `parameters.hzeeman=(Z+Z')/2; clear('ZI','ZQ')`.
- Lines 123: computes `C` using `C=carrier(spin_system,parameters.rframes{k}{1})`.
- Lines 132: computes `R` using `R=relaxation(spin_system,parameters.orientation)`.
- Lines 135: computes `K` using `K=kinetics(spin_system)`.
- Lines 138: computes `parameters.spc_dim` using `parameters.spc_dim=1; parameters.spn_dim=size(H,1)`.
- Lines 144: computes `answer` using `answer=pulse_sequence(spin_system,parameters,H,R,K)`.

### Local helper functions

- Line 149: `defaults()` — `function parameters=defaults(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'decouple')`.
  - Representative operation: `report(spin_system,'parameters.decouple field not set, assuming no decoupling.')`.
- Line 173: `grumble()` — `function grumble(spin_system,pulse_sequence,parameters,assumptions)`. Orientation
  - Representative operation: `if ~isfield(parameters,'orientation')`.
  - Representative operation: `error('system orientation must be specified in parameters.orientation variable.')`.

## Parameters / inputs

- pulse_sequence -a function handle to one of the pulse se-
- quences located in the experiments folder
- assumptions -a string that would be passed to assume.m
- when the Hamiltonian is built
- parameters.spins -a cell array giving the spins that
- the pulse sequence involves, e.g.
- {'1H','13C'}
- parameters.offset -a cell array giving transmitter off-
- sets in Hz on each of the spins listed
- in parameters.spins array
- parameters.orientation -a row vector of the three Euler angles
- (in radians) giving the orientation of
- the system relative to the input orien-
- tation.
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N,3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.needs -a cell array of strings specifying additional
- information required by the sequence:
- 'zeeman_op' -Zeeman part of the Hamiltonian
- in the laboratory frame, to be placed into
- parameters.hzeeman and sent to pulse sequence
- 'aniso_eq' -thermal equilibrium is recomputed
- using the full anisotropic Hamiltonian at the
- current orientation, and sent to the pulse
- sequence in parameters.rho0 subfield
- parameters.* -additional subfields may be required by your
- pulse sequence -check its documentation page
- The parameters structure is passed to the pulse sequence with the follo-
- wing additional parameters set:
- parameters.spc_dim -matrix dimension for the spatial
- dynamics subspace (1 in this case)
- parameters.spn_dim -matrix dimension for the spin
- dynamics subspace

## Outputs

- this function returns whatever it is that the pulse sequence returns
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.

## Implementation structure

- Single-crystal interface to pulse sequences. Generates a Liouvillian
- superoperator and passes it on to the pulse sequence function, which
- should be supplied as a handle. Syntax:
- answer=crystal(spin_system,pulse_sequence,parameters,assumptions)
- pulse_sequence -a function handle to one of the pulse se-
- quences located in the experiments folder
- assumptions -a string that would be passed to assume.m
- when the Hamiltonian is built
- parameters.spins -a cell array giving the spins that
- the pulse sequence involves, e.g.
- {'1H','13C'}
- parameters.offset -a cell array giving transmitter off-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `report()`, `assume()`, `hamiltonian()`, `ismember()`, `equilibrium()`, `frqoffset()`, `orientation()`, `clear()`, `carrier()`, `rotframe()`, `relaxation()`, `kinetics()`, `pulse_sequence()`.
