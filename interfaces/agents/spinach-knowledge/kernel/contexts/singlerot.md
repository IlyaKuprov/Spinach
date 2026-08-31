# kernel/contexts/singlerot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/singlerot.m`
- Signature: `[answer,sph_grid]=singlerot(spin_system,pulse_sequence,...`
- Total lines: 522

## Purpose

Single angle spinning context. In Liouville space, this wrapper builds the Fokker-Planck evolution generator that includes the spin Hamiltoni- an commutation superoperator, applicable dissipators (relaxation, kine- tics), and the rotor turning generator. In Hilbert space, this wrapper builds the stack of spin Hamiltonians, one for each rotor phase. Those are handed over to the pulse sequence, which the user must supp

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `parfor_progr()`, `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 104-105: Show the banner; implemented by `banner(spin_system,'sequence_banner')`.
- Lines 107-108: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 110-111: Check consistency; implemented by `grumble(spin_system,pulse_sequence,parameters,assumptions)`.
- Lines 113-114: Set the interaction assumptions; implemented by `spin_system=assume(spin_system,assumptions)`.
- Lines 116-117: Get isotropic Hamiltonian and rotations basis; implemented by `[I,Q]=hamiltonian(assume(spin_system,assumptions))`.
- Lines 119-120: Apply channel frequency offsets; implemented by `I=frqoffset(spin_system,I,parameters)`.
- Lines 122-123: Compute isotropic thermal equilibrium; implemented by `if ismember('iso_eq',parameters.needs)`.
- Lines 129-131: Get carrier operators for numerical rotating frame transformations; implemented by `C=cell(size(parameters.rframes))`.
- Lines 136-137: Get relaxation and kinetics generators; implemented by `R=relaxation(spin_system); K=kinetics(spin_system)`.
- Lines 139-141: Load the spherical integration grid; implemented by `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas','gammas','weights')`.
- Lines 145-148: Get and report rotor axis orientation angles; implemented by `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 153-154: Report problem dimensions to the user; implemented by `spn_dim=sum(size(I,2))`.
- Lines 160-161: Forward dimension information to the pulse sequence; implemented by `parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim`.
- Lines 163-164: Formalism-dependent stage; implemented by `switch spin_system.bas.formalism`.
- Lines 166-167: Liouville space; implemented by `case {'sphten-liouv','zeeman-liouv'}`.
- Lines 169-170: Report the composite dimension of the Fokker-Planck problem; implemented by `report(spin_system,['Fokker-Planck problem dimension: ' num2str(spc_dim*spn_dim)])`.
- Lines 172-173: Make the rotor turning generator; implemented by `[rotor_phases,d_dphi]=fourdif(spc_dim,1)`.
- Lines 176-177: Project relaxation and kinetics superoperators into the FP space; implemented by `R=kron(speye([spc_dim spc_dim]),R); K=kron(speye([spc_dim spc_dim]),K)`.

### Control flow inferred from the code

- Line 123: conditional branch on `ismember('iso_eq',parameters.needs)`.
- Line 132: `for` loop over `n=1:numel(parameters.rframes)`.
- Line 164: dispatches on `spin_system.bas.formalism`; cases `{'sphten-liouv','zeeman-liouv'}`, `{'zeeman-hilb','zeeman-wavef'}`.
- Line 180: conditional branch on `isfield(parameters,'rho0')&&strcmp(parameters.grid,'single_crystal')`.
- Line 197: conditional branch on `isfield(parameters,'coil')`.
- Line 219: conditional branch on `isfield(parameters,'serial')&&parameters.serial`.
- Line 233: conditional branch on `~isfield(parameters,'verbose')||(parameters.verbose==0)`.
- Line 239: conditional branch on `~isworkernode`.

### Key state/data transformations

- Lines 108: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 114: computes `spin_system` using `spin_system=assume(spin_system,assumptions)`.
- Lines 117: computes `[I,Q]` using `[I,Q]=hamiltonian(assume(spin_system,assumptions))`.
- Lines 120: computes `I` using `I=frqoffset(spin_system,I,parameters)`.
- Lines 125: computes `I_labframe` using `I_labframe=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 126: computes `parameters.rho0` using `parameters.rho0=equilibrium(spin_system,I_labframe)`.
- Lines 131: computes `C` using `C=cell(size(parameters.rframes))`.
- Lines 133: computes `C{n}` using `C{n}=carrier(spin_system,parameters.rframes{n}{1})`.
- Lines 137: computes `R` using `R=relaxation(spin_system); K=kinetics(spin_system)`.
- Lines 140-141: computes `sph_grid` using `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas','gammas','weights')`.
- Lines 142: computes `alphas` using `alphas=sph_grid.alphas; betas=sph_grid.betas`.
- Lines 143: computes `gammas` using `gammas=sph_grid.gammas; weights=sph_grid.weights`.
- Lines 146-148: computes `[rotor_phi,rotor_theta,~]` using `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 149: computes `rotor_theta` using `rotor_theta=pi/2-rotor_theta`.
- Lines 154: computes `spn_dim` using `spn_dim=sum(size(I,2))`.
- Lines 155: computes `spc_dim` using `spc_dim=2*parameters.max_rank+1; n_orients=numel(weights)`.
- Lines 161: computes `parameters.spc_dim` using `parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim`.
- Lines 173: computes `[rotor_phases,d_dphi]` using `[rotor_phases,d_dphi]=fourdif(spc_dim,1)`.

### Local helper functions

- Line 249: `parfor_progr()` — `function parfor_progr()`. Parallel powder averaging loop
  - Representative operation: `orients_done=orients_done+1; last_message=toc-last_toc`.
  - Representative operation: `if (last_message>5)||(orients_done==n_orients)`.
- Line 375: `defaults()` — `function parameters=defaults(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'decouple')`.
  - Representative operation: `report(spin_system,'parameters.decouple field not set, assuming no decoupling.')`.
- Line 401: `grumble()` — `function grumble(spin_system,pulse_sequence,parameters,assumptions)`. Option combination restrictions
  - Representative operation: `if (~iscell(parameters.needs))||any(~cellfun(@ischar,parameters.needs(:)))`.
  - Representative operation: `error('parameters.needs must be a cell array of character strings.')`.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.rate -spinning rate in Hz. Positive numbers
- for JEOL, negative for Varian and Bruker
- due to different rotation directions.
- parameters.axis -spinning axis, given as a normalized
- 3-element vector; this is the direction
- around which the rotor is turning
- parameters.spins -a cell array giving the spins that
- the pulse sequence involves, e.g.
- {'1H','13C'}
- parameters.offset -a cell array giving transmitter off-
- sets in Hz on each of the spins listed
- in parameters.spins array
- parameters.max_rank -maximum rotor harmonic rank to retain
- in the solution (increase till conver-
- gence is achieved, a good guess value
- is the number of spinning sidebands
- expected in the spectrum)
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N',3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.grid -spherical grid file name; see grids
- directory in the kernel. Two-angle
- grids should be used in Liouville
- space and three-angle grids in Hil-
- bert space.
- parameters.needs -a cell array of character strings spe-
- cifying additional requirements that
- the sequence has:
- 'iso_eq' -thermal equilibrium state
- of the isotropic Hamiltonian will be
- placed into parameters.rho0
- parameters.sum_up -when set to 1 (default), returns the
- powder average. When set to 0, returns
- individual answers for each point in
- the powder as a cell array.
- parameters.* -additional subfields may be required
- by the pulse sequence -check its do-
- cumentation page
- assumptions -context-specific assumptions ('nmr', 'epr',
- 'labframe', etc.) -see the pulse sequence
- header for information on this setting.
- The parameters structure is passed to the pulse sequence with the follo-
- wing additional parameters set:
- parameters.spc_dim -matrix dimension for the spatial
- dynamics subspace
- parameters.spn_dim -matrix dimension for the spin
- dynamics subspace

## Outputs

- answer -the poweder average or a cell array ofwhatever it is
- that the pulse sequence returns
- sph_grid -spherical grid used ithe calculation
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.

## Implementation structure

- Single angle spinning context. In Liouville space, this wrapper builds
- the Fokker-Planck evolution generator that includes the spin Hamiltoni-
- an commutation superoperator, applicable dissipators (relaxation, kine-
- tics), and the rotor turning generator. In Hilbert space, this wrapper
- builds the stack of spin Hamiltonians, one for each rotor phase. Those
- are handed over to the pulse sequence, which the user must supply as a
- function handle. Syntax:
- [answer,sph_grid]=singlerot(spin_system,pulse_sequence,...
- parameters,assumptions)
- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `assume()`, `hamiltonian()`, `frqoffset()`, `ismember()`, `report()`, `equilibrium()`, `carrier()`, `relaxation()`, `kinetics()`, `load()`, `cart2sph()`, `num2str()`, `fourdif()`.
