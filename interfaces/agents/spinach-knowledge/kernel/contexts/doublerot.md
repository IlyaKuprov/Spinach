# kernel/contexts/doublerot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/doublerot.m`
- Signature: `[answer,sph_grid]=doublerot(spin_system,pulse_sequence,...`
- Total lines: 551

## Purpose

Double angle spinning context. In Liouville space, this wrapper builds the Fokker-Planck evolution generator and passes it on to the pulse se- quence function, which should be supplied as a handle. In Hilbert space, this wrapper builds the stack of spin Hamiltonians, one for each pair of rotor phases on the two-rotor phase grid, and hands that stack to the pulse sequence. Syntax: [answer,sph_grid]=doublerot(spin_syst

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 111-112: Show the banner; implemented by `banner(spin_system,'sequence_banner')`.
- Lines 114-115: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 117-118: Check consistency; implemented by `grumble(spin_system,pulse_sequence,parameters,assumptions)`.
- Lines 120-121: Report to the user; implemented by `report(spin_system,'building the evolution generator ')`.
- Lines 123-124: Set the assumptions; implemented by `spin_system=assume(spin_system,assumptions)`.
- Lines 126-127: Get the Hamiltonian; implemented by `[H,Q]=hamiltonian(spin_system)`.
- Lines 129-130: Apply offsets; implemented by `H=frqoffset(spin_system,H,parameters)`.
- Lines 132-133: Compute the thermal equilibrium; implemented by `if ismember('iso_eq',parameters.needs)`.
- Lines 141-142: Count rotor trajectory points; implemented by `npoints_inner=2*parameters.rank_inner+1`.
- Lines 146-147: Get problem dimensions; implemented by `spc_dim=npoints_total; spn_dim=size(H,1)`.
- Lines 152-153: Compute spectral derivative operators; implemented by `[traj_inner,d_dphi_inner]=fourdif(npoints_inner,1)`.
- Lines 156-157: Compute rotor phase tracks; implemented by `phases_outer=kron(traj_outer,ones(npoints_inner,1))`.
- Lines 160-161: Formalism-dependent stage; implemented by `switch spin_system.bas.formalism`.
- Lines 163-164: Liouville space; implemented by `case {'sphten-liouv','zeeman-liouv'}`.
- Lines 166-167: Report the composite dimension of the Fokker-Planck problem; implemented by `report(spin_system,['Fokker-Planck problem dimension ' num2str(spc_dim*spn_dim)])`.
- Lines 169-171: Compute double spinning operator; implemented by `M=2*pi*kron((parameters.rate_outer*kron(d_dphi_outer,speye(size(d_dphi_inner)))+ parameters.rate_inner*kron(speye(size(d_dphi_outer)),d_dphi_inner)),speye(size(H)))`.
- Lines 173-174: Hilbert space; implemented by `case {'zeeman-hilb','zeeman-wavef'}`.
- Lines 176-177: No spinning operator in the rotor stack route; implemented by `M=[]`.

### Control flow inferred from the code

- Line 133: conditional branch on `ismember('iso_eq',parameters.needs)`.
- Line 135: conditional branch on `isfield(parameters,'rho0')`.
- Line 161: dispatches on `spin_system.bas.formalism`; cases `{'sphten-liouv','zeeman-liouv'}`, `{'zeeman-hilb','zeeman-wavef'}`.
- Line 197: `for` loop over `n=1:numel(parameters.rframes)`.
- Line 213: conditional branch on `ismember(spin_system.bas.formalism,{'zeeman-hilb','zeeman-wavef'})&&`.
- Line 219: conditional branch on `ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
- Line 225: conditional branch on `isfield(parameters,'rho0')&&strcmp(parameters.grid,'single_crystal')`.
- Line 242: conditional branch on `isfield(parameters,'coil')`.
- Line 253: conditional branch on `isfield(parameters,'serial')&&`.
- Line 272: conditional branch on `~isfield(parameters,'verbose')||(parameters.verbose==0)`.
- Line 281: `parfor` loop over `(q=1:numel(weights),nworkers)`.
- Line 285: `for` loop over `n=1:npoints_total`.
- Line 286: `for` loop over `k=1:npoints_total`.
- Line 292: `for` loop over `n=1:npoints_total`.

### Key state/data transformations

- Lines 115: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 124: computes `spin_system` using `spin_system=assume(spin_system,assumptions)`.
- Lines 127: computes `[H,Q]` using `[H,Q]=hamiltonian(spin_system)`.
- Lines 130: computes `H` using `H=frqoffset(spin_system,H,parameters)`.
- Lines 138: computes `parameters.rho0` using `parameters.rho0=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 142: computes `npoints_inner` using `npoints_inner=2*parameters.rank_inner+1`.
- Lines 143: computes `npoints_outer` using `npoints_outer=2*parameters.rank_outer+1`.
- Lines 144: computes `npoints_total` using `npoints_total=npoints_outer*npoints_inner`.
- Lines 147: computes `spc_dim` using `spc_dim=npoints_total; spn_dim=size(H,1)`.
- Lines 150: computes `parameters.spc_dim` using `parameters.spc_dim=spc_dim; parameters.spn_dim=spn_dim`.
- Lines 153: computes `[traj_inner,d_dphi_inner]` using `[traj_inner,d_dphi_inner]=fourdif(npoints_inner,1)`.
- Lines 154: computes `[traj_outer,d_dphi_outer]` using `[traj_outer,d_dphi_outer]=fourdif(npoints_outer,1)`.
- Lines 157: computes `phases_outer` using `phases_outer=kron(traj_outer,ones(npoints_inner,1))`.
- Lines 158: computes `phases_inner` using `phases_inner=kron(ones(npoints_outer,1),traj_inner)`.
- Lines 170-171: computes `M` using `M=2*pi*kron((parameters.rate_outer*kron(d_dphi_outer,speye(size(d_dphi_inner)))+ parameters.rate_inner*kron(speye(size(d_dphi_outer)),d_dphi_inner)),speye(size(H)))`.
- Lines 187-189: computes `[phi_inner,theta_inner,~]` using `[phi_inner,theta_inner,~]=cart2sph(parameters.axis_inner(1), parameters.axis_inner(2), parameters.axis_inner(3))`.
- Lines 190-192: computes `[phi_outer,theta_outer,~]` using `[phi_outer,theta_outer,~]=cart2sph(parameters.axis_outer(1), parameters.axis_outer(2), parameters.axis_outer(3))`.
- Lines 193: computes `theta_inner` using `theta_inner=pi/2-theta_inner; theta_outer=pi/2-theta_outer`.

### Local helper functions

- Line 398: `defaults()` — `function parameters=defaults(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'decouple')`.
  - Representative operation: `report(spin_system,'parameters.decouple field not set, assuming no decoupling.')`.
- Line 424: `grumble()` — `function grumble(spin_system,pulse_sequence,parameters,assumptions)`. Wavefunctions cannot represent thermal equilibria
  - Representative operation: `if strcmp(spin_system.bas.formalism,'zeeman-wavef')&&isfield(parameters,'needs')&& ismember('iso_eq',parameters.needs)`.
  - Representative operation: `ismember('iso_eq',parameters.needs)`.

## Outputs

- answer -the poweder average or a cell array ofwhatever it is
- that the pulse sequence returns
- sph_grid -spherical grid used ithe calculation
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.
- Note: the state projector assumes a powder --single crystal DOR is not
- currently supported.
- Note: the function supports parallel processing via Matlab's Distri-
- buted Computing Toolbox -different system orientations are eva-
- luated on different labs.

## Implementation structure

- Double angle spinning context. In Liouville space, this wrapper builds
- the Fokker-Planck evolution generator and passes it on to the pulse se-
- quence function, which should be supplied as a handle. In Hilbert space,
- this wrapper builds the stack of spin Hamiltonians, one for each pair of
- rotor phases on the two-rotor phase grid, and hands that stack to the
- pulse sequence. Syntax:
- [answer,sph_grid]=doublerot(spin_system,pulse_sequence,...
- parameters,assumptions)
- where pulse sequence is a function handle to one of the pulse sequences
- located in the experiments directory, assumptions is a string that would
- be passed to assume.m when the Hamiltonian is built and parameters is a
- structure with the following subfields:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `report()`, `assume()`, `hamiltonian()`, `frqoffset()`, `ismember()`, `isfield()`, `equilibrium()`, `num2str()`, `fourdif()`, `speye()`, `cart2sph()`, `carrier()`, `relaxation()`.
