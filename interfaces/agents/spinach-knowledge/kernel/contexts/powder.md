# kernel/contexts/powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/powder.m`
- Signature: `[answer,sph_grid]=powder(spin_system,pulse_sequence,...`
- Total lines: 441

## Purpose

Static powder interface to pulse sequences. Generates a Liouvillian superoperator, the initial state, the coil state, then passes them on to the pulse sequence function. Syntax: [answer,sph_grid]=powder(spin_system,pulse_sequence,... parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `parfor_progr()`, `defaults()`, `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 102-103: Show the banner; implemented by `banner(spin_system,'sequence_banner')`.
- Lines 105-106: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 108-109: Check consistency; implemented by `grumble(spin_system,pulse_sequence,parameters,assumptions)`.
- Lines 111-112: Report to the user; implemented by `report(spin_system,'building the Liouvillian ')`.
- Lines 114-115: Get the lab frame Zeeman operator if needed; implemented by `if ismember('zeeman_op',parameters.needs)`.
- Lines 117-118: Isotropic + anisotropic basis, for assembly later; implemented by `report(spin_system,'building the lab frame Zeeman operator ')`.
- Lines 123-124: Work around the parfor bug; implemented by `ZI=[]; ZQ=[]`.
- Lines 128-129: Get the lab frame Hamiltonian if needed; implemented by `if ismember('iso_eq',parameters.needs)`.
- Lines 131-132: Isotropic part of the lab frame Hamiltonian; implemented by `report(spin_system,'building the lab frame Hamiltonian ')`.
- Lines 135-136: Compute thermal equilibrium here; implemented by `parameters.rho0=equilibrium(spin_system,HL); QL=[]`.
- Lines 140-141: Isotropic + anisotropic part, for thermal equilibrium later; implemented by `report(spin_system,'building the lab frame Hamiltonian ')`.
- Lines 146-147: Work around the parfor bug; implemented by `HL=[]; QL=[]`.
- Lines 151-152: Set the assumptions specified by the user; implemented by `spin_system=assume(spin_system,assumptions)`.
- Lines 154-155: Get Hamiltonian building blocks; implemented by `[I,Q]=hamiltonian(spin_system)`.
- Lines 157-158: Get kinetics superoperator; implemented by `K=kinetics(spin_system)`.
- Lines 160-161: Add offsets to the isotropic part; implemented by `I=frqoffset(spin_system,I,parameters)`.
- Lines 163-164: Get carrier operators; implemented by `C=cell(size(parameters.rframes))`.
- Lines 169-170: Get problem dimensions; implemented by `parameters.spc_dim=1; parameters.spn_dim=size(I,1)`.

### Control flow inferred from the code

- Line 115: conditional branch on `ismember('zeeman_op',parameters.needs)`.
- Line 129: conditional branch on `ismember('iso_eq',parameters.needs)`.
- Line 165: `for` loop over `n=1:numel(parameters.rframes)`.
- Line 185: conditional branch on `isfield(parameters,'serial')&&parameters.serial`.
- Line 200: conditional branch on `~isfield(parameters,'verbose')||(parameters.verbose==0)`.
- Line 206: conditional branch on `~isworkernode`.

### Key state/data transformations

- Lines 106: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 119: computes `[ZI,ZQ]` using `[ZI,ZQ]=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 124: computes `ZI` using `ZI=[]; ZQ=[]`.
- Lines 133: computes `HL` using `HL=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 136: computes `parameters.rho0` using `parameters.rho0=equilibrium(spin_system,HL); QL=[]`.
- Lines 142: computes `[HL,QL]` using `[HL,QL]=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 152: computes `spin_system` using `spin_system=assume(spin_system,assumptions)`.
- Lines 155: computes `[I,Q]` using `[I,Q]=hamiltonian(spin_system)`.
- Lines 158: computes `K` using `K=kinetics(spin_system)`.
- Lines 161: computes `I` using `I=frqoffset(spin_system,I,parameters)`.
- Lines 164: computes `C` using `C=cell(size(parameters.rframes))`.
- Lines 166: computes `C{n}` using `C{n}=carrier(spin_system,parameters.rframes{n}{1})`.
- Lines 170: computes `parameters.spc_dim` using `parameters.spc_dim=1; parameters.spn_dim=size(I,1)`.
- Lines 173-174: computes `sph_grid` using `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas','gammas','weights')`.
- Lines 177: computes `alphas` using `alphas=sph_grid.alphas; betas=sph_grid.betas`.
- Lines 178: computes `gammas` using `gammas=sph_grid.gammas; weights=sph_grid.weights`.
- Lines 179: computes `n_orients` using `n_orients=numel(weights)`.
- Lines 182: computes `ans_array` using `ans_array=cell(n_orients,1)`.

### Local helper functions

- Line 216: `parfor_progr()` — `function parfor_progr()`. Parallel powder averaging loop
  - Representative operation: `orients_done=orients_done+1; last_message=toc-last_toc`.
  - Representative operation: `if (last_message>5)||(orients_done==n_orients)`.
- Line 312: `defaults()` — `function parameters=defaults(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'decouple')`.
  - Representative operation: `report(spin_system,'parameters.decouple field not set, assuming no decoupling.')`.
- Line 340: `grumble()` — `function grumble(spin_system,pulse_sequence,parameters,assumptions)`. Spherical grid
  - Representative operation: `if ~isfield(parameters,'grid')`.
  - Representative operation: `error('spherical averaging grid must be specified in parameters.grid variable.')`.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the spins that the
- pulse sequence works on, in the order
- of channels, e.g. {'1H','13C'}
- parameters.offset -a cell array giving transmitter offsets
- in Hz on each of the spins listed in
- parameters.spins
- parameters.grid -name of the spherical averaging grid
- file (see the grids directory in the
- kernel).
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N,3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.needs -a cell array of strings specifying ad-
- ditional information required by the
- sequence:
- 'zeeman_op' -Zeeman part of the Hami-
- ltonian in the laboratory frame, to be
- placed into parameters.hzeeman and sent
- to the pulse sequence
- 'iso_eq' -thermal equilibrium is com-
- computed using the isotropic part of
- the Hamiltonian, and sent to the pulse
- sequence via parameters.rho0
- 'aniso_eq' -thermal equilibrium is re-
- computed using the full anisotropic Ha-
- miltonian at each orientation, and sent
- to pulse sequence via parameters.rho0
- parameters.rho0 -initial state; may be a function handle
- that depends on the three Euler angles
- in ZYZ active convention
- parameters.serial -if set to true, disables automatic pa-
- rallelisation
- parameters.sum_up -if set to false, causes the pulse sequ-
- ence output at each orientation to be
- returned instead of the powder average
- parameters.* -additional subfields may be required by
- the pulse sequence -check its documen-
- tation page
- assumptions -context-specific assumptions ('nmr', 'epr',
- 'labframe', etc.) -see the pulse sequence
- header for information on this setting.

## Outputs

- answer -powder average of whatever it is that the pulse
- sequence returns; if parameters.sum_up is set to
- false, a cell array of outputs at each orienta-
- tion is returned
- sph_grid -powder averaging grid data structure with three
- Euler angles and weights for each point
- Note: THIS IS FOR STATIC POWDERS -use singlerot for MAS simulations.
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.
- Note: the function supports parallel processing via Matlab's Distri-
- buted Computing Toolbox -different system orientations are eva-
- luated on different labs.

## Implementation structure

- Static powder interface to pulse sequences. Generates a Liouvillian
- superoperator, the initial state, the coil state, then passes them
- on to the pulse sequence function. Syntax:
- [answer,sph_grid]=powder(spin_system,pulse_sequence,...
- parameters,assumptions)
- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the spins that the
- pulse sequence works on, in the order
- of channels, e.g. {'1H','13C'}
- parameters.offset -a cell array giving transmitter offsets

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `report()`, `ismember()`, `hamiltonian()`, `assume()`, `equilibrium()`, `kinetics()`, `frqoffset()`, `carrier()`, `load()`, `isfield()`, `num2str()`, `afterEach()`, `ticBytes()`.
