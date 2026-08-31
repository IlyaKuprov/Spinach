# kernel/contexts/floquet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/floquet.m`
- Signature: `[answer,sph_grid]=floquet(spin_system,pulse_sequence,...`
- Total lines: 433

## Purpose

Floquet magic angle spinning context. Generates a Liouvillian super- operator and passes it on to the pulse sequence function, which sho- uld be supplied as a handle. Syntax: [answer,sph_grid]=floquet(spin_system,pulse_sequence,... parameters,assumptions) where pulse sequence is a function handle to one of the pulse sequences located in the experiments directory, assumptions is a string that would be passed to assume

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `parfor_progr()`, `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 86-87: Show the banner; implemented by `banner(spin_system,'sequence_banner')`.
- Lines 89-90: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 92-93: Check consistency; implemented by `grumble(spin_system,pulse_sequence,parameters,assumptions)`.
- Lines 95-96: Report to the user; implemented by `report(spin_system,'building the Liouvillian ')`.
- Lines 98-99: Set the assumptions; implemented by `spin_system=assume(spin_system,assumptions)`.
- Lines 101-102: Get the Hamiltonian; implemented by `[H,Q]=hamiltonian(spin_system)`.
- Lines 104-105: Detect non-empty spherical ranks; implemented by `int_ranks=[]`.
- Lines 113-114: Fourier block extent covers every non-empty rank; implemented by `max_int=max([int_ranks 2])`.
- Lines 116-117: Warn about truncated rotor harmonics; implemented by `if max_int>parameters.max_rank`.
- Lines 121-122: Apply frequency offsets; implemented by `H=frqoffset(spin_system,H,parameters)`.
- Lines 124-125: Compute the thermal equilibrium; implemented by `if ismember('iso_eq',parameters.needs)`.
- Lines 133-134: Get problem dimensions; implemented by `spc_dim=2*parameters.max_rank+1; spn_dim=size(H,1)`.
- Lines 140-141: Get relaxation and kinetics; implemented by `R=relaxation(spin_system); K=kinetics(spin_system)`.
- Lines 143-145: Get the averaging grid as a structure; implemented by `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas','gammas','weights')`.
- Lines 147-148: Assign local variables; implemented by `alphas=sph_grid.alphas; betas=sph_grid.betas`.
- Lines 152-153: Rotor turning generator, spinning sense matched to singlerot.m; implemented by `M=2*pi*parameters.rate*kron(spdiags((-parameters.max_rank:parameters.max_rank)',0,spc_dim,spc_dim),speye(spn_dim))`.
- Lines 155-156: Kron up relaxation and kinetics; implemented by `R=kron(speye(spc_dim),R); K=kron(speye(spc_dim),K)`.
- Lines 158-159: Get the rotor orientation; implemented by `[phi,theta,~]=cart2sph(parameters.axis(1),parameters.axis(2),parameters.axis(3))`.

### Control flow inferred from the code

- Line 106: `for` loop over `r=1:numel(Q)`.
- Line 107: conditional branch on `any(cellfun(@nnz,Q{r}),'all')`.
- Line 117: conditional branch on `max_int>parameters.max_rank`.
- Line 125: conditional branch on `ismember('iso_eq',parameters.needs)`.
- Line 127: conditional branch on `isfield(parameters,'rho0')`.
- Line 164: `for` loop over `r=int_ranks`.
- Line 170: conditional branch on `isfield(parameters,'rho0')`.
- Line 173: conditional branch on `isfield(parameters,'coil')`.
- Line 184: conditional branch on `isfield(parameters,'serial')&&`.
- Line 203: conditional branch on `~isfield(parameters,'verbose')||(parameters.verbose==0)`.
- Line 209: conditional branch on `~isworkernode`.

### Key state/data transformations

- Lines 90: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 99: computes `spin_system` using `spin_system=assume(spin_system,assumptions)`.
- Lines 102: computes `[H,Q]` using `[H,Q]=hamiltonian(spin_system)`.
- Lines 105: computes `int_ranks` using `int_ranks=[]`.
- Lines 114: computes `max_int` using `max_int=max([int_ranks 2])`.
- Lines 122: computes `H` using `H=frqoffset(spin_system,H,parameters)`.
- Lines 130: computes `parameters.rho0` using `parameters.rho0=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 134: computes `spc_dim` using `spc_dim=2*parameters.max_rank+1; spn_dim=size(H,1)`.
- Lines 138: computes `parameters.spc_dim` using `parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim`.
- Lines 141: computes `R` using `R=relaxation(spin_system); K=kinetics(spin_system)`.
- Lines 144-145: computes `sph_grid` using `sph_grid=load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'alphas','betas','gammas','weights')`.
- Lines 148: computes `alphas` using `alphas=sph_grid.alphas; betas=sph_grid.betas`.
- Lines 149: computes `gammas` using `gammas=sph_grid.gammas; weights=sph_grid.weights`.
- Lines 150: computes `n_orients` using `n_orients=numel(weights)`.
- Lines 153: computes `M` using `M=2*pi*parameters.rate*kron(spdiags((-parameters.max_rank:parameters.max_rank)',0,spc_dim,spc_dim),speye(spn_dim))`.
- Lines 159: computes `[phi,theta,~]` using `[phi,theta,~]=cart2sph(parameters.axis(1),parameters.axis(2),parameters.axis(3))`.
- Lines 160: computes `theta` using `theta=pi/2-theta`.
- Lines 163: computes `D_lab2rot` using `D_lab2rot=cell(1,max_int); D_rot2lab=cell(1,max_int)`.

### Local helper functions

- Line 219: `parfor_progr()` — `function parfor_progr()`. Powder averaged spectrum
  - Representative operation: `orients_done=orients_done+1; last_message=toc-last_toc`.
  - Representative operation: `if (last_message>5)||(orients_done==n_orients)`.
- Line 318: `defaults()` — `function parameters=defaults(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'decouple')`.
  - Representative operation: `report(spin_system,'parameters.decouple field not set, assuming no decoupling.')`.
- Line 340: `grumble()` — `function grumble(spin_system,pulse_sequence,parameters,assumptions)`. Rotating frames
  - Representative operation: `if isfield(parameters,'rframes')`.
  - Representative operation: `error('numerical rotating frame transformation is not supported by Floquet theory.')`.

## Outputs

- answer -the poweder average or a cell array ofwhatever it is
- that the pulse sequence returns
- sph_grid -spherical grid used ithe calculation
- Note: the choice of the rank depends on the spinning rate (the slower
- the spinning, the greater ranks are required). The rank is appro-
- ximately equal to the number of spinning sidebands.
- Note: the state projector assumes a powder --single crystal MAS is not
- currently supported.
- Note: perturbative corrections to the rotating frame transformation are
- not supported -use singlerot.m instead.
- Note: the spinning sense matches singlerot.m -the same parameters.rate
- produces the same powder result in both contexts.
- Note: the function supports parallel processing via Matlab's Distri-
- buted Computing Toolbox -different system orientations are eva-
- luated on different labs.

## Implementation structure

- Floquet magic angle spinning context. Generates a Liouvillian super-
- operator and passes it on to the pulse sequence function, which sho-
- uld be supplied as a handle. Syntax:
- [answer,sph_grid]=floquet(spin_system,pulse_sequence,...
- parameters,assumptions)
- where pulse sequence is a function handle to one of the pulse sequences
- located in the experiments directory, assumptions is a string that would
- be passed to assume.m when the Hamiltonian is built and parameters is a
- structure with the following subfields:
- parameters.rate -spinning rate in Hz. Positive numbers
- for JEOL, negative for Varian and Bruker
- due to different rotation directions.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `report()`, `assume()`, `hamiltonian()`, `any()`, `cellfun()`, `num2str()`, `frqoffset()`, `ismember()`, `isfield()`, `equilibrium()`, `relaxation()`, `kinetics()`, `load()`.
