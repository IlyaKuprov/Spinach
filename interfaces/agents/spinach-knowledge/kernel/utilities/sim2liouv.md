# kernel/utilities/sim2liouv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sim2liouv.m`
- Signature: `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`
- Total lines: 159

## Purpose

Moves a zeeman-hilb simulation context into Liouville space. When the formalism specified in the spin system object is 'zeeman-hilb', this function projects the evolution generators into Liouville space, converts the standard state-like and operator-like fields of the parameters structure, rebuilds the basis index table, mig- rates the symmetry irrep projectors into the adjoint representa- tion, and sets the formalis

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 64-65: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 67-68: Only the zeeman-hilb formalism needs the move; implemented by `if strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Lines 70-71: Inform the user; implemented by `report(spin_system,'projecting zeeman-hilb simulation into Liouville space ')`.
- Lines 73-74: Project the evolution generators into Liouville space; implemented by `H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm')`.
- Lines 76-77: Get the Hilbert space basis table and dimension; implemented by `zbas=spin_system.bas.basis; hdim=size(zbas,1)`.
- Lines 79-80: Stretch the state-like parameters; implemented by `if isfield(parameters,'rho0')`.
- Lines 90-91: Project the operator-like parameters; implemented by `if isfield(parameters,'pulse_op')`.
- Lines 101-102: Rebuild the basis index table for the Liouville space; implemented by `spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))]`.
- Lines 104-105: Migrate the irreps into the adjoint representation; implemented by `if isfield(spin_system.bas,'irrep')`.
- Lines 107-108: Grab the Hilbert space irreps; implemented by `hs_irreps=spin_system.bas.irrep; n_irreps=numel(hs_irreps)`.
- Lines 110-111: Preallocate the Liouville space irrep array; implemented by `ls_irreps(n_irreps^2)=struct('projector',[],'dimension',[])`.
- Lines 113-114: Loop over ordered pairs of Hilbert space irreps; implemented by `for n=1:n_irreps`.
- Lines 117-118: Build the irrep pair projector and dimension; implemented by `pair_idx=n_irreps*(n-1)+k`.
- Lines 127-128: Write the Liouville space irreps; implemented by `spin_system.bas.irrep=ls_irreps`.
- Lines 134-135: Update the formalism setting; implemented by `spin_system.bas.formalism='zeeman-liouv'`.

### Control flow inferred from the code

- Line 68: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Line 80: conditional branch on `isfield(parameters,'rho0')`.
- Line 83: conditional branch on `isfield(parameters,'coil')`.
- Line 86: conditional branch on `isfield(parameters,'screen')`.
- Line 91: conditional branch on `isfield(parameters,'pulse_op')`.
- Line 94: conditional branch on `isfield(parameters,'mw_oper')`.
- Line 97: conditional branch on `isfield(parameters,'ez_oper')`.
- Line 105: conditional branch on `isfield(spin_system.bas,'irrep')`.
- Line 114: `for` loop over `n=1:n_irreps`.
- Line 115: `for` loop over `k=1:n_irreps`.

### Key state/data transformations

- Lines 74: computes `H` using `H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm')`.
- Lines 77: computes `zbas` using `zbas=spin_system.bas.basis; hdim=size(zbas,1)`.
- Lines 81: computes `parameters.rho0` using `parameters.rho0=reshape(parameters.rho0,hdim^2,[])`.
- Lines 84: computes `parameters.coil` using `parameters.coil=reshape(parameters.coil,hdim^2,[])`.
- Lines 87: computes `parameters.screen` using `parameters.screen=reshape(parameters.screen,hdim^2,[])`.
- Lines 92: computes `parameters.pulse_op` using `parameters.pulse_op=hilb2liouv(parameters.pulse_op,'comm')`.
- Lines 95: computes `parameters.mw_oper` using `parameters.mw_oper=hilb2liouv(parameters.mw_oper,'comm')`.
- Lines 98: computes `parameters.ez_oper` using `parameters.ez_oper=hilb2liouv(parameters.ez_oper,'comm')`.
- Lines 102: computes `spin_system.bas.basis` using `spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))]`.
- Lines 108: computes `hs_irreps` using `hs_irreps=spin_system.bas.irrep; n_irreps=numel(hs_irreps)`.
- Lines 111: computes `ls_irreps(n_irreps^2)` using `ls_irreps(n_irreps^2)=struct('projector',[],'dimension',[])`.
- Lines 118: computes `pair_idx` using `pair_idx=n_irreps*(n-1)+k`.
- Lines 119-120: computes `ls_irreps(pair_idx).projector` using `ls_irreps(pair_idx).projector=kron(conj(hs_irreps(n).projector), hs_irreps(k).projector)`.
- Lines 121-122: computes `ls_irreps(pair_idx).dimension` using `ls_irreps(pair_idx).dimension=hs_irreps(n).dimension* hs_irreps(k).dimension`.
- Lines 128: computes `spin_system.bas.irrep` using `spin_system.bas.irrep=ls_irreps`.
- Lines 135: computes `spin_system.bas.formalism` using `spin_system.bas.formalism='zeeman-liouv'`.

### Local helper functions

- Line 142: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv', 'zeeman-hilb','zeeman-wavef'})`.
  - Representative operation: `'zeeman-hilb','zeeman-wavef'})`.

## Parameters / inputs

- spin_system -Spinach spin system object
- parameters -pulse sequence parameters structure; the
- state-like fields rho0, coil, and screen
- (matrices or their horizontal concatena-
- tions) are stretched into state vectors,
- and the operator-like fields pulse_op,
- mw_oper, and ez_oper are converted into
- commutation superoperators, when present
- H -Hamiltonian operator, converted into a
- commutation superoperator; an empty
- matrix is passed through
- R -relaxation matrix, converted into an
- anticommutation superoperator; an empty
- matrix is passed through
- K -kinetics matrix, converted into an
- anticommutation superoperator; an empty
- matrix is passed through

## Outputs

- spin_system -spin system object with zeeman-liouv
- formalism and basis information
- parameters -parameters structure with the standard
- fields converted into Liouville space
- H,R,K -Liouville space evolution generators
- Note: a Hilbert space density matrix block S(n)*Y*S(k)' maps to
- the state vector kron(conj(S(k)),S(n))*Y(:), and so each
- ordered pair of Hilbert space irrep projectors yields the
- Liouville space irrep projector kron(conj(S(k)),S(n)).
- Every such subspace is invariant under superoperators
- built from symmetry-respecting Hilbert space generators;
- unpopulated subspaces are dropped by reduce.m at run time
- in the usual way.

## Implementation structure

- Moves a zeeman-hilb simulation context into Liouville space. When
- the formalism specified in the spin system object is 'zeeman-hilb',
- this function projects the evolution generators into Liouville
- space, converts the standard state-like and operator-like fields
- of the parameters structure, rebuilds the basis index table, mig-
- rates the symmetry irrep projectors into the adjoint representa-
- tion, and sets the formalism to 'zeeman-liouv'; for all other
- formalisms, every argument is returned unchanged. This makes
- Liouville-space pulse sequences callable with zeeman-hilb
- inputs. Syntax:
- [spin_system,parameters,H,R,K]=...
- sim2liouv(spin_system,parameters,H,R,K)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `report()`, `hilb2liouv()`, `isfield()`, `ls_irreps()`, `conj()`, `hs_irreps()`, `num2str()`, `ismember()`, `isstruct()`.
