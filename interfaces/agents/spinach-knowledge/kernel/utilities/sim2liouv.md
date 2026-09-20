# kernel/utilities/sim2liouv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sim2liouv.m`
- Signature: `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`
- Total lines: 190

## Purpose

Moves a zeeman-hilb simulation context into Liouville space. When the formalism specified in the spin system object is 'zeeman-hilb', this function projects the evolution generators into Liouville space, converts the standard state-like and operator-like fields of the parameters structure, rebuilds the basis index table, mig- rates the symmetry irrep projectors into the adjoint representa- tion, and sets the formalis

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The anticommutation superoperator of a Hilbert space relaxation matrix damps the unit state, which the Liouville space branch of relaxation.m never does. The row and the column of the unit state are therefore projected out of the converted R, so that the unit state is neither damped nor a source of relaxation and the trace is conserved. For the scalar damping matrix that the kernel builds in Hilbert space the result coincides with the Liouville space damp operator. The unit state spans all diagonal irrep pairs, and those are merged into one subspace so that the projected R stays block-diagonal in the irrep table that reduce.m evolves independently.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.
- The unit state exemption is the symmetric projection `R=R-U*(U'*R)-(R*U)*U'+U*(U'*R*U)*U'` with `U=unit_state(spin_system)` taken after the formalism switch, i.e. the normalised stretched unit matrix; it is exact for any relaxation matrix, scalar or not, and leaves R Hermitian when R is Hermitian.
- With symmetry, the migrated irrep table has `n_irreps^2-n_irreps+1` entries: the first holds all diagonal irrep pairs `kron(conj(S(n)),S(n))` side by side (dimension is the sum of the squared irrep dimensions), the rest are the off-diagonal pairs. The unit state lives entirely in the first entry, so the projection does not leak between reduction blocks, and population contrasts between irreps are damped exactly as in the native zeeman-liouv damp operator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 78-79: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 81-82: Only the zeeman-hilb formalism needs the move; implemented by `if strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Lines 84-85: Inform the user; implemented by `report(spin_system,'projecting zeeman-hilb simulation into Liouville space ')`.
- Lines 87-88: Project the evolution generators into Liouville space; implemented by `H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm')`.
- Lines 90-91: Get the Hilbert space basis table and dimension; implemented by `zbas=spin_system.bas.basis; hdim=size(zbas,1)`.
- Lines 93-94: Stretch the state-like parameters; implemented by `if isfield(parameters,'rho0')`.
- Lines 104-105: Project the operator-like parameters; implemented by `if isfield(parameters,'pulse_op')`.
- Lines 115-116: Rebuild the basis index table for the Liouville space; implemented by `spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))]`.
- Lines 118-119: Migrate the irreps into the adjoint representation; implemented by `if isfield(spin_system.bas,'irrep')`.
- Lines 121-122: Grab the Hilbert space irreps; implemented by `hs_irreps=spin_system.bas.irrep; n_irreps=numel(hs_irreps)`.
- Lines 124-125: Preallocate the Liouville space irrep array; implemented by `ls_irreps(n_irreps^2-n_irreps+1)=struct('projector',[],'dimension',[])`.
- Lines 127-128: Diagonal irrep pairs share the unit state and are merged into the first subspace; implemented by `ls_irreps(1).projector=[]; ls_irreps(1).dimension=0; pair_idx=1`.
- Lines 130-131: Loop over ordered pairs of Hilbert space irreps; implemented by `for n=1:n_irreps`.
- Lines 134-136: Build the irrep pair projector and dimension; implemented by `pair_proj=kron(conj(hs_irreps(n).projector),hs_irreps(k).projector); pair_dim=hs_irreps(n).dimension*hs_irreps(k).dimension`.
- Lines 138-139: Merge diagonal pairs, store off-diagonal pairs separately; implemented by `if n==k`.
- Lines 151-152: Write the Liouville space irreps; implemented by `spin_system.bas.irrep=ls_irreps`.
- Lines 158-159: Update the formalism setting; implemented by `spin_system.bas.formalism='zeeman-liouv'`.
- Lines 161-165: Project the unit state out of the relaxation superoperator; implemented by `U=unit_state(spin_system); R=R-U*(U'*R)-(R*U)*U'+U*(U'*R*U)*U'` followed by `report(spin_system,'unit state exempted from the projected relaxation superoperator.')`.

### Control flow inferred from the code

- Line 82: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Line 94: conditional branch on `isfield(parameters,'rho0')`.
- Line 97: conditional branch on `isfield(parameters,'coil')`.
- Line 100: conditional branch on `isfield(parameters,'screen')`.
- Line 105: conditional branch on `isfield(parameters,'pulse_op')`.
- Line 108: conditional branch on `isfield(parameters,'mw_oper')`.
- Line 111: conditional branch on `isfield(parameters,'ez_oper')`.
- Line 119: conditional branch on `isfield(spin_system.bas,'irrep')`.
- Line 131: `for` loop over `n=1:n_irreps`.
- Line 132: `for` loop over `k=1:n_irreps`.
- Line 139: conditional branch on `n==k`, merging diagonal pairs into `ls_irreps(1)` and storing off-diagonal pairs at `pair_idx`.
- Line 162: conditional branch on `~isempty(R)`; an empty relaxation matrix is passed through unchanged.

### Key state/data transformations

- Lines 88: computes `H` using `H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm')`.
- Lines 91: computes `zbas` using `zbas=spin_system.bas.basis; hdim=size(zbas,1)`.
- Lines 95: computes `parameters.rho0` using `parameters.rho0=reshape(parameters.rho0,hdim^2,[])`.
- Lines 98: computes `parameters.coil` using `parameters.coil=reshape(parameters.coil,hdim^2,[])`.
- Lines 101: computes `parameters.screen` using `parameters.screen=reshape(parameters.screen,hdim^2,[])`.
- Lines 106: computes `parameters.pulse_op` using `parameters.pulse_op=hilb2liouv(parameters.pulse_op,'comm')`.
- Lines 109: computes `parameters.mw_oper` using `parameters.mw_oper=hilb2liouv(parameters.mw_oper,'comm')`.
- Lines 112: computes `parameters.ez_oper` using `parameters.ez_oper=hilb2liouv(parameters.ez_oper,'comm')`.
- Lines 116: computes `spin_system.bas.basis` using `spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))]`.
- Lines 122: computes `hs_irreps` using `hs_irreps=spin_system.bas.irrep; n_irreps=numel(hs_irreps)`.
- Lines 125: computes `ls_irreps(n_irreps^2-n_irreps+1)` using `ls_irreps(n_irreps^2-n_irreps+1)=struct('projector',[],'dimension',[])`.
- Lines 128: computes `ls_irreps(1)` and `pair_idx` using `ls_irreps(1).projector=[]; ls_irreps(1).dimension=0; pair_idx=1`.
- Lines 135-136: computes `pair_proj` and `pair_dim` using `pair_proj=kron(conj(hs_irreps(n).projector),hs_irreps(k).projector); pair_dim=hs_irreps(n).dimension*hs_irreps(k).dimension`.
- Lines 140-141: computes `ls_irreps(1).projector` and `ls_irreps(1).dimension` using `ls_irreps(1).projector=[ls_irreps(1).projector pair_proj]; ls_irreps(1).dimension=ls_irreps(1).dimension+pair_dim`.
- Lines 143-145: computes `pair_idx`, `ls_irreps(pair_idx).projector`, and `ls_irreps(pair_idx).dimension` using `pair_idx=pair_idx+1; ls_irreps(pair_idx).projector=pair_proj; ls_irreps(pair_idx).dimension=pair_dim`.
- Lines 152: computes `spin_system.bas.irrep` using `spin_system.bas.irrep=ls_irreps`.
- Lines 159: computes `spin_system.bas.formalism` using `spin_system.bas.formalism='zeeman-liouv'`.
- Lines 163-164: computes `U` and `R` using `U=unit_state(spin_system); R=R-U*(U'*R)-(R*U)*U'+U*(U'*R*U)*U'`.

### Local helper functions

- Line 173: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
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
- anticommutation superoperator with the
- unit state exempted from damping; an
- empty matrix is passed through
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
- Note: the anticommutation superoperator of a Hilbert space rela-
- xation matrix damps the unit state, which the Liouville
- space branch of relaxation.m never does. The row and the
- column of the unit state are therefore projected out of
- the converted R, so that the unit state is neither damped
- nor a source of relaxation and the trace is conserved. For
- the scalar damping matrix that the kernel builds in Hilbert
- space the result coincides with the Liouville space damp
- operator. The unit state spans all diagonal irrep pairs,
- and those are merged into one subspace so that the projec-
- ted R stays block-diagonal in the irrep table that reduce.m
- evolves independently.

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

- Called routines detected from the main body: `grumble()`, `strcmp()`, `report()`, `hilb2liouv()`, `isfield()`, `ls_irreps()`, `conj()`, `hs_irreps()`, `numel()`, `num2str()`, `unit_state()`, `isempty()`, `ismember()`, `isstruct()`.
