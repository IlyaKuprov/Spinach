# kernel/contexts/gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/contexts/gridfree.m`
- Signature: `answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)`
- Total lines: 365

## Purpose

Fokker-Planck magic angle spinning and SLE context. Generates a Liouvil- lian superoperator and passes it on to the pulse sequence function, which should be supplied as a handle. Syntax: answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 81-82: Show the banner; implemented by `banner(spin_system,'sequence_banner')`.
- Lines 84-85: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 87-88: Check consistency; implemented by `grumble(spin_system,pulse_sequence,parameters,assumptions)`.
- Lines 90-91: Report to the user; implemented by `report(spin_system,'building the Liouvillian ')`.
- Lines 93-94: Set the assumptions; implemented by `spin_system=assume(spin_system,assumptions)`.
- Lines 96-97: Get the Hamiltonian; implemented by `[H,Q]=hamiltonian(spin_system)`.
- Lines 99-100: Detect non-empty spherical ranks; implemented by `int_ranks=[]`.
- Lines 108-109: Warn about spatially truncated interaction ranks; implemented by `if any(int_ranks>parameters.max_rank)`.
- Lines 113-114: Get spatial operators; implemented by `[Lx,Ly,Lz,D,~]=sle_operators(parameters.max_rank,int_ranks)`.
- Lines 116-117: Apply offsets; implemented by `H=frqoffset(spin_system,H,parameters)`.
- Lines 119-120: Add user-specified operators; implemented by `if isfield(parameters,'add_terms')`.
- Lines 127-128: Compute isotropic thermal equilibrium; implemented by `if ismember('iso_eq',parameters.needs)`.
- Lines 136-137: Get relaxation and kinetics; implemented by `R=relaxation(spin_system)`.
- Lines 140-141: Get problem dimensions; implemented by `spc_dim=size(Lx,1); parameters.spc_dim=spc_dim`.
- Lines 144-145: Inform the user; implemented by `report(spin_system,['spin space problem dimension ' num2str(spn_dim)])`.
- Lines 149-150: Project isotropic parts; implemented by `H=polyadic({{speye(spc_dim),H}})`.
- Lines 154-155: Add anisotropic parts of every non-empty rank; implemented by `for r=int_ranks`.
- Lines 163-164: Rotation operator; implemented by `if isfield(parameters,'rate')&&abs(parameters.rate)>0`.

### Control flow inferred from the code

- Line 101: `for` loop over `r=1:numel(Q)`.
- Line 102: conditional branch on `any(cellfun(@nnz,Q{r}),'all')`.
- Line 109: conditional branch on `any(int_ranks>parameters.max_rank)`.
- Line 120: conditional branch on `isfield(parameters,'add_terms')`.
- Line 121: `for` loop over `n=1:numel(parameters.add_terms)`.
- Line 128: conditional branch on `ismember('iso_eq',parameters.needs)`.
- Line 130: conditional branch on `isfield(parameters,'rho0')`.
- Line 155: `for` loop over `r=int_ranks`.
- Line 156: `for` loop over `k=1:(2*r+1)`.
- Line 157: `for` loop over `m=1:(2*r+1)`.
- Line 164: conditional branch on `isfield(parameters,'rate')&&abs(parameters.rate)>0`.
- Line 182: conditional branch on `isfield(parameters,'tau_c')`.
- Line 185: conditional branch on `isscalar(parameters.tau_c)`.
- Line 213: conditional branch on `~ismember('polyadic',spin_system.sys.enable)`.

### Key state/data transformations

- Lines 85: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 94: computes `spin_system` using `spin_system=assume(spin_system,assumptions)`.
- Lines 97: computes `[H,Q]` using `[H,Q]=hamiltonian(spin_system)`.
- Lines 100: computes `int_ranks` using `int_ranks=[]`.
- Lines 114: computes `[Lx,Ly,Lz,D,~]` using `[Lx,Ly,Lz,D,~]=sle_operators(parameters.max_rank,int_ranks)`.
- Lines 117: computes `H` using `H=frqoffset(spin_system,H,parameters)`.
- Lines 133: computes `parameters.rho0` using `parameters.rho0=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 137: computes `R` using `R=relaxation(spin_system)`.
- Lines 138: computes `K` using `K=kinetics(spin_system)`.
- Lines 141: computes `spc_dim` using `spc_dim=size(Lx,1); parameters.spc_dim=spc_dim`.
- Lines 142: computes `spn_dim` using `spn_dim=size(H,1); parameters.spn_dim=spn_dim`.
- Lines 167: computes `spinning_axis` using `spinning_axis=parameters.axis/norm(parameters.axis,2)`.
- Lines 174: computes `Hr` using `Hr=2*pi*parameters.rate*(spinning_axis(1)*Lx+spinning_axis(2)*Ly+spinning_axis(3)*Lz)`.
- Lines 188: computes `diff_iso` using `diff_iso=1/(6*parameters.tau_c)`.
- Lines 189: computes `Rd` using `Rd=diff_iso*(Lx*Lx+Ly*Ly+Lz*Lz)`.
- Lines 195: computes `dten` using `dten=inv(6*parameters.tau_c)`.
- Lines 220: computes `H_whos` using `H_whos=whos('H'); report(spin_system,['memory footprint of H array: ' num2str(H_whos.bytes/1024^3) ' GB'])`.
- Lines 221: computes `R_whos` using `R_whos=whos('R'); report(spin_system,['memory footprint of R array: ' num2str(R_whos.bytes/1024^3) ' GB'])`.

### Local helper functions

- Line 242: `defaults()` — `function parameters=defaults(spin_system,parameters)`.
  - Representative operation: `if ~isfield(parameters,'offset')`.
  - Representative operation: `report(spin_system,'parameters.offset field not set, assuming zero offsets.')`.
- Line 257: `grumble()` — `function grumble(spin_system,pulse_sequence,parameters,assumptions)`. Rotating frames
  - Representative operation: `if isfield(parameters,'rframes')`.
  - Representative operation: `error('numerical rotating frame transformation is not supported by SLE formalism.')`.

## Parameters / inputs

- pulse_sequence -a function handle to one of the pulse sequences
- located in the experiments directory
- assumptions -is a string that would be passed to assume.m
- when the Hamiltonian is built
- parameters -a structure with the following subfields:
- .rate -spinning rate in Hz. Positive numbers
- for JEOL, negative for Varian and Bruker
- due to different rotation directions.
- .axis -spinning axis, given as a normalized
- 3-element vector
- .spins -a cell array giving the spins that
- the pulse sequence involves, e.g.
- {'1H','13C'}
- .offset -a cell array giving transmitter off-
- sets in Hz on each of the spins listed
- in parameters.spins array
- .max_rank -maximum D-function rank to retain in
- the solution (increase till conver-
- gence is achieved, approximately
- equal to the number of spinning si-
- debands in the spectrum)
- .tau_c -correlation times (in seconds) for rotational
- diffusion. Single number for isotropic rotati-
- onal diffusion, and a symmetric positive defi-
- nite 3x3 correlation time tensor for anisotro-
- pic rotational diffusion; the rotational dif-
- fusion tensor is inv(6*tau_c).
- .* -additional subfields may be required by your
- pulse sequence -check its documentation page
- The parameters structure is passed to the pulse sequence with the follo-
- wing additional parameters set:
- parameters.spc_dim -matrix dimension for the spatial
- dynamics subspace
- parameters.spn_dim -matrix dimension for the spin
- dynamics subspace

## Outputs

- this context function returns the powder average of whatever it
- is that the pulse sequence returns
- Note: the choice of the Wigner D function rank truncation level depends on
- on the spinning rate (the slower the spinning, the greater ranks are
- required).
- Note: rotational correlation times for SLE go into parameters.tau_c, not
- inter.tau_c (the latter is only used by the Redfield theory module).
- Note: the state projector assumes a powder --single crystal MAS is not
- currently supported.
- Note: perturbative corrections to the rotating frame transformation are
- not supported -use singlerot.m if you need them.

## Implementation structure

- Fokker-Planck magic angle spinning and SLE context. Generates a Liouvil-
- lian superoperator and passes it on to the pulse sequence function, which
- should be supplied as a handle. Syntax:
- answer=gridfree(spin_system,pulse_sequence,parameters,assumptions)
- pulse_sequence -a function handle to one of the pulse sequences
- located in the experiments directory
- assumptions -is a string that would be passed to assume.m
- when the Hamiltonian is built
- parameters -a structure with the following subfields:
- .rate -spinning rate in Hz. Positive numbers
- for JEOL, negative for Varian and Bruker
- due to different rotation directions.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `banner()`, `defaults()`, `grumble()`, `report()`, `assume()`, `hamiltonian()`, `any()`, `cellfun()`, `num2str()`, `sle_operators()`, `frqoffset()`, `isfield()`, `ismember()`, `equilibrium()`, `relaxation()`, `kinetics()`.
