# kernel/utilities/sim2liouv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sim2liouv.m`
- Signature: `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`
- Total lines: 194

## Purpose

Moves a zeeman-hilb simulation context into Liouville space. When the formalism specified in the spin system object is 'zeeman-hilb', this function projects the evolution generators into Liouville space, converts the standard state-like and operator-like fields of the parameters structure, rebuilds the basis index table, mig- rates the symmetry irrep projectors into the adjoint representa- tion, and sets the formalis

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The anticommutation superoperator of a Hilbert space relaxation matrix damps the unit state, which the Liouville space branch of relaxation.m never does. The row and the column of the unit state are therefore projected out of the converted R, so that the unit state is neither damped nor a source of relaxation and the trace is conserved. With symmetry, this is done for the unit state of every irrep, which keeps R block-diagonal in the irrep pair subspaces that reduce.m evolves independently. For the scalar damping matrix that the kernel builds in Hilbert space the result coincides with the Liouville space damp operator.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.
- The unit state exemption is the symmetric projection `R=R-U*(U'*R)-(R*U)*U'+U*(U'*R*U)*U'`, where the columns of `U` are the normalised stretched projectors of the Hilbert space irreps (or the normalised stretched unit matrix without symmetry); the columns are orthonormal, so `U*U'` is the orthogonal projector onto the block unit states, and the result is exact for any relaxation matrix, scalar or not.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 77-78: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 80-81: Only the zeeman-hilb formalism needs the move; implemented by `if strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Lines 83-84: Inform the user; implemented by `report(spin_system,'projecting zeeman-hilb simulation into Liouville space ')`.
- Lines 86-87: Project the evolution generators into Liouville space; implemented by `H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm')`.
- Lines 89-90: Get the Hilbert space basis table and dimension; implemented by `zbas=spin_system.bas.basis; hdim=size(zbas,1)`.
- Lines 92-93: Stretch the state-like parameters; implemented by `if isfield(parameters,'rho0')`.
- Lines 103-104: Project the operator-like parameters; implemented by `if isfield(parameters,'pulse_op')`.
- Lines 114-115: Rebuild the basis index table for the Liouville space; implemented by `spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))]`.
- Lines 117-118: Migrate the irreps into the adjoint representation; implemented by `if isfield(spin_system.bas,'irrep')`.
- Lines 120-121: Grab the Hilbert space irreps; implemented by `hs_irreps=spin_system.bas.irrep; n_irreps=numel(hs_irreps)`.
- Lines 123-124: Preallocate the Liouville space irrep array; implemented by `ls_irreps(n_irreps^2)=struct('projector',[],'dimension',[])`.
- Lines 126-127: Preallocate the irrep unit state array; implemented by `unit_cols=cell(1,n_irreps)`.
- Lines 129-130: Loop over ordered pairs of Hilbert space irreps; implemented by `for n=1:n_irreps`.
- Lines 132-134: Normalised unit state of the irrep; implemented by `unit_n=sparse(reshape(hs_irreps(n).projector*hs_irreps(n).projector',[],1)); unit_cols{n}=unit_n/norm(unit_n,2)`.
- Lines 138-139: Build the irrep pair projector and dimension; implemented by `pair_idx=n_irreps*(n-1)+k`.
- Lines 148-149: Assemble the irrep unit state array; implemented by `U=[unit_cols{:}]`.
- Lines 151-152: Write the Liouville space irreps; implemented by `spin_system.bas.irrep=ls_irreps`.
- Lines 158-159: Normalised unit state of the whole space; implemented by `U=reshape(speye(hdim),[],1)/sqrt(hdim)`.
- Lines 163-164: Update the formalism setting; implemented by `spin_system.bas.formalism='zeeman-liouv'`.
- Lines 166-169: Project the block unit states out of the relaxation superoperator; implemented by `R=R-U*(U'*R)-(R*U)*U'+U*(U'*R*U)*U'` followed by `report(spin_system,'unit state exempted from the projected relaxation superoperator.')`.

### Control flow inferred from the code

- Line 81: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Line 93: conditional branch on `isfield(parameters,'rho0')`.
- Line 96: conditional branch on `isfield(parameters,'coil')`.
- Line 99: conditional branch on `isfield(parameters,'screen')`.
- Line 104: conditional branch on `isfield(parameters,'pulse_op')`.
- Line 107: conditional branch on `isfield(parameters,'mw_oper')`.
- Line 110: conditional branch on `isfield(parameters,'ez_oper')`.
- Line 118: conditional branch on `isfield(spin_system.bas,'irrep')`, with an `else` branch at line 156 for systems without symmetry.
- Line 130: `for` loop over `n=1:n_irreps`.
- Line 136: `for` loop over `k=1:n_irreps`.
- Line 167: conditional branch on `~isempty(R)`; an empty relaxation matrix is passed through unchanged.

### Key state/data transformations

- Lines 87: computes `H` using `H=hilb2liouv(H,'comm'); R=hilb2liouv(R,'acomm'); K=hilb2liouv(K,'acomm')`.
- Lines 90: computes `zbas` using `zbas=spin_system.bas.basis; hdim=size(zbas,1)`.
- Lines 94: computes `parameters.rho0` using `parameters.rho0=reshape(parameters.rho0,hdim^2,[])`.
- Lines 97: computes `parameters.coil` using `parameters.coil=reshape(parameters.coil,hdim^2,[])`.
- Lines 100: computes `parameters.screen` using `parameters.screen=reshape(parameters.screen,hdim^2,[])`.
- Lines 105: computes `parameters.pulse_op` using `parameters.pulse_op=hilb2liouv(parameters.pulse_op,'comm')`.
- Lines 108: computes `parameters.mw_oper` using `parameters.mw_oper=hilb2liouv(parameters.mw_oper,'comm')`.
- Lines 111: computes `parameters.ez_oper` using `parameters.ez_oper=hilb2liouv(parameters.ez_oper,'comm')`.
- Lines 115: computes `spin_system.bas.basis` using `spin_system.bas.basis=[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))]`.
- Lines 121: computes `hs_irreps` using `hs_irreps=spin_system.bas.irrep; n_irreps=numel(hs_irreps)`.
- Lines 124: computes `ls_irreps(n_irreps^2)` using `ls_irreps(n_irreps^2)=struct('projector',[],'dimension',[])`.
- Lines 127: computes `unit_cols` using `unit_cols=cell(1,n_irreps)`.
- Lines 133-134: computes `unit_cols{n}` using `unit_n=sparse(reshape(hs_irreps(n).projector*hs_irreps(n).projector',[],1)); unit_cols{n}=unit_n/norm(unit_n,2)`.
- Lines 139: computes `pair_idx` using `pair_idx=n_irreps*(n-1)+k`.
- Lines 140-141: computes `ls_irreps(pair_idx).projector` using `ls_irreps(pair_idx).projector=kron(conj(hs_irreps(n).projector), hs_irreps(k).projector)`.
- Lines 142-143: computes `ls_irreps(pair_idx).dimension` using `ls_irreps(pair_idx).dimension=hs_irreps(n).dimension* hs_irreps(k).dimension`.
- Lines 149: computes `U` using `U=[unit_cols{:}]`.
- Lines 152: computes `spin_system.bas.irrep` using `spin_system.bas.irrep=ls_irreps`.
- Lines 159: computes `U` using `U=reshape(speye(hdim),[],1)/sqrt(hdim)`.
- Lines 164: computes `spin_system.bas.formalism` using `spin_system.bas.formalism='zeeman-liouv'`.
- Lines 168: computes `R` using `R=R-U*(U'*R)-(R*U)*U'+U*(U'*R*U)*U'`.

### Local helper functions

- Line 177: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
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
- nor a source of relaxation and the trace is conserved. With
- symmetry, this is done for the unit state of every irrep,
- which keeps R block-diagonal in the irrep pair subspaces
- that reduce.m evolves independently. For the scalar dam-
- ping matrix that the kernel builds in Hilbert space the
- result coincides with the Liouville space damp operator.

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

- Called routines detected from the main body: `grumble()`, `strcmp()`, `report()`, `hilb2liouv()`, `isfield()`, `ls_irreps()`, `conj()`, `hs_irreps()`, `num2str()`, `sparse()`, `reshape()`, `norm()`, `speye()`, `sqrt()`, `isempty()`, `ismember()`, `isstruct()`.
