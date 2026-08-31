# kernel/contexts/device.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/device.m`
- Signature: `answer=device(spin_system,pulse_sequence,parameters,assumptions)`
- Total lines: 276

## Purpose

Spin-boson device interface to pulse sequences. Generates the evolution generators for a device containing spins and bosonic modes at a fixed orientation of the spin subsystem, and passes them to the pulse sequence function, which should be supplied as a handle. Syntax: answer=device(spin_system,pulse_sequence,parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 74-75: Show the banner; implemented by `banner(spin_system,'sequence_banner')`.
- Lines 77-78: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 80-81: Check consistency; implemented by `grumble(spin_system,pulse_sequence,parameters,assumptions)`.
- Lines 83-84: Set assumptions; implemented by `spin_system=assume(spin_system,assumptions)`.
- Lines 86-87: Get the Hamiltonian at the spin subsystem orientation; implemented by `[I,Q]=hamiltonian(spin_system)`.
- Lines 90-91: Get relaxation and kinetics superoperators; implemented by `R=relaxation(spin_system,parameters.orientation)`.
- Lines 94-95: Get the thermal equilibrium if needed; implemented by `if ismember('rho_eq',parameters.needs)`.
- Lines 101-102: Process spin channel offsets; implemented by `if ~isempty(parameters.spins)`.
- Lines 106-107: Process mode detuning offsets; implemented by `mode_list=find(ismember(spin_system.comp.types,{'C','V','T'}))`.
- Lines 116-117: Get carrier operators; implemented by `C=cell(size(parameters.rframes))`.
- Lines 122-123: Apply rotating frames; implemented by `for k=1:numel(parameters.rframes)`.
- Lines 127-128: Wipe the states of the decoupled spins; implemented by `H=decouple(spin_system,H,[],parameters.decouple)`.
- Lines 135-136: Get problem dimensions; implemented by `parameters.spc_dim=1; parameters.spn_dim=size(H,1)`.
- Lines 138-139: Report to the user; implemented by `report(spin_system,'running the pulse sequence ')`.
- Lines 141-142: Call the pulse sequence; implemented by `answer=pulse_sequence(spin_system,parameters,H,R,K)`.

### Control flow inferred from the code

- Line 95: conditional branch on `ismember('rho_eq',parameters.needs)`.
- Line 102: conditional branch on `~isempty(parameters.spins)`.
- Line 108: `for` loop over `n=1:numel(mode_list)`.
- Line 109: conditional branch on `abs(parameters.mode_offset(n))>0`.
- Line 118: `for` loop over `n=1:numel(parameters.rframes)`.
- Line 123: `for` loop over `k=1:numel(parameters.rframes)`.
- Line 131: conditional branch on `ismember('rho_eq',parameters.needs)`.

### Key state/data transformations

- Lines 78: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 84: computes `spin_system` using `spin_system=assume(spin_system,assumptions)`.
- Lines 87: computes `[I,Q]` using `[I,Q]=hamiltonian(spin_system)`.
- Lines 88: computes `H` using `H=I+orientation(Q,parameters.orientation)`.
- Lines 91: computes `R` using `R=relaxation(spin_system,parameters.orientation)`.
- Lines 92: computes `K` using `K=kinetics(spin_system)`.
- Lines 97: computes `[HL,QL]` using `[HL,QL]=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 98: computes `parameters.rho0` using `parameters.rho0=equilibrium(spin_system,HL,QL,parameters.orientation)`.
- Lines 107: computes `mode_list` using `mode_list=find(ismember(spin_system.comp.types,{'C','V','T'}))`.
- Lines 117: computes `C` using `C=cell(size(parameters.rframes))`.
- Lines 119: computes `C{n}` using `C{n}=carrier(spin_system,parameters.rframes{n}{1})`.
- Lines 132: computes `[~,parameters.rho0]` using `[~,parameters.rho0]=decouple(spin_system,[],parameters.rho0,parameters.decouple)`.
- Lines 136: computes `parameters.spc_dim` using `parameters.spc_dim=1; parameters.spn_dim=size(H,1)`.
- Lines 142: computes `answer` using `answer=pulse_sequence(spin_system,parameters,H,R,K)`.

### Local helper functions

- Line 147: `defaults()` — `function parameters=defaults(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'spins')`.
  - Representative operation: `report(spin_system,'parameters.spins field not set, assuming no spin channels.')`.
- Line 179: `grumble()` — `function grumble(spin_system,pulse_sequence,parameters,assumptions)`.
  - Representative operation: `if ~isa(pulse_sequence,'function_handle')`.
  - Representative operation: `error('pulse_sequence argument must be a function handle.')`.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the spin species that
- the pulse sequence works on, in the order
- of channels, e.g. {'E'}; may be omitted
- when no spin channels are needed
- parameters.offset -transmitter offsets in Hz on each of the
- spin species listed in parameters.spins
- parameters.mode_offset -detuning offsets in Hz, one for each
- bosonic mode in the order of declaration;
- the transmitter sign convention of para-
- meters.offset applies: each offset enters
- the Hamiltonian as minus the offset times
- the number operator of its mode, so that a
- positive offset lowers the mode frequency
- parameters.decouple -a cell array of spin species to be wiped
- from the evolution generators and from the
- initial state, e.g. {'1H'}; the default is
- an empty cell array, meaning no decoupling
- parameters.orientation -Euler angles (ZYZ active convention, ra-
- dians) giving the orientation of the spin
- subsystem; bosonic terms are not affected;
- the default is [0 0 0]
- parameters.rframes -numerical rotating frame specification for
- spin species, e.g. {{'E',2}}, as described
- in the header of rotframe.m
- parameters.needs -a cell array of strings specifying additional
- information required by the sequence:
- 'rho_eq' -thermal equilibrium state at the
- system temperature, with the Bose-Einstein
- populations of the bosonic modes included;
- this is placed into parameters.rho0
- parameters.* -additional subfields may be required by
- the pulse sequence -check its docs
- assumptions -'labframe', 'cavity', or 'spin-phonon';
- see the header of assume.m for details

## Outputs

- answer -whatever it is that the pulse sequence returns.
- Note: the system must contain at least one bosonic mode; pure spin
- systems belong with liquid.m, crystal.m, and powder.m contexts.
- Note: dissipative bosonic modes require a Liouville space formalism;
- coherent simulations may also use zeeman-hilb.

## Implementation structure

- Spin-boson device interface to pulse sequences. Generates the evolution
- generators for a device containing spins and bosonic modes at a fixed
- orientation of the spin subsystem, and passes them to the pulse sequence
- function, which should be supplied as a handle. Syntax:
- answer=device(spin_system,pulse_sequence,parameters,assumptions)
- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the spin species that
- the pulse sequence works on, in the order
- of channels, e.g. {'E'}; may be omitted
- when no spin channels are needed

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `assume()`, `hamiltonian()`, `orientation()`, `relaxation()`, `kinetics()`, `ismember()`, `report()`, `equilibrium()`, `frqoffset()`, `num2str()`, `mode_list()`, `operator()`, `carrier()`.
