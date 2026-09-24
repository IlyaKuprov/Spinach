# kernel/utilities/sim2liouv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/sim2liouv.m`
- Signature: `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`
- Total lines: 193

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

## Parameters / inputs

- spin_system -Spinach spin system object
- parameters -pulse sequence parameters structure; the
- state-like fields rho0, coil, and screen
- (matrices or their horizontal concatena-
- tions) are stretched into state vectors,
- and the operator-like fields pulse_op,
- mw_oper, ez_oper, and homodec_oper are converted into
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
