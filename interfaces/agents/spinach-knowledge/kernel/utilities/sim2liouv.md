# kernel/utilities/sim2liouv.m

- Signature: `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`

## Purpose

Moves a zeeman-hilb simulation context into Liouville space. When the spin system formalism is 'zeeman-hilb', this function projects the evolution generators into Liouville space, converts the standard state-like and operator-like fields of the parameters structure, rebuilds the basis index table, migrates the symmetry irrep projectors into the adjoint representation, and sets the formalism to 'zeeman-liouv'. Other formalisms return every argument unchanged, allowing Liouville-space pulse sequences to be called with zeeman-hilb inputs.

## Physical / mathematical content

- The anticommutation superoperator of a Hilbert space relaxation matrix damps the unit state, unlike the Liouville-space branch of relaxation.m. Projecting the unit-state row and column out of the converted R makes the unit state neither damped nor a source of relaxation and conserves the trace. For the scalar damping matrix built by the kernel in Hilbert space, the result coincides with the native zeeman-liouv damp operator. The unit state spans all diagonal irrep pairs, which are merged into one subspace so that the projected R remains block-diagonal in the irrep table evolved independently by reduce.m.
- A Hilbert-space density matrix block `S(n)*Y*S(k)'` maps to `kron(conj(S(k)),S(n))*Y(:)`; each ordered pair of Hilbert-space irrep projectors therefore yields a Liouville-space projector `kron(conj(S(k)),S(n))`. These subspaces are invariant under superoperators built from symmetry-respecting Hilbert-space generators. reduce.m drops unpopulated subspaces at run time.

## Numerical / algorithmic content

- The converted basis table refreshes existing cache metadata using the canonical `md5_hash` rule, keeping Hilbert operators and Hamiltonians separate from their Liouville representations. Objects without cache metadata gain none; non-Hilbert inputs remain unchanged.
- The unit-state exemption uses the symmetric projection `R=R-U*(U'*R)-(R*U)*U'+U*(U'*R*U)*U'`, where `U=unit_state(spin_system)` is taken after the formalism switch and is the normalised stretched unit matrix. For any relaxation matrix, scalar or not, the projection makes both the unit-state row and column exactly zero and preserves Hermiticity when R is Hermitian.
- With symmetry, the migrated irrep table has `n_irreps^2-n_irreps+1` entries. The first combines all diagonal irrep-pair projectors `kron(conj(S(n)),S(n))` side by side, with dimension equal to the sum of their squared dimensions; the remaining entries hold off-diagonal pairs. The unit state lies entirely in the first entry, so the projection does not leak between reduction blocks. Population contrasts between irreps are damped exactly as in the native zeeman-liouv damp operator.

## Parameters / inputs

- spin_system - Spinach spin system object.
- parameters - Pulse sequence parameters structure. When present, state-like fields rho0, coil, and screen (matrices or their horizontal concatenations) are stretched into state vectors; operator-like fields pulse_op, mw_oper, ez_oper, and homodec_oper become commutation superoperators.
- H - Hamiltonian operator converted into a commutation superoperator; an empty matrix passes through.
- R - Relaxation matrix converted into an anticommutation superoperator with the unit state exempted from damping; an empty matrix passes through.
- K - Kinetics matrix converted into an anticommutation superoperator; an empty matrix passes through.

## Outputs

- spin_system - Spin system object with zeeman-liouv formalism and basis information.
- parameters - Parameters structure with the standard fields converted into Liouville space.
- H,R,K - Liouville-space evolution generators.

## Header notes

Only zeeman-hilb inputs are converted; other formalisms return every input unchanged.

ilya.kuprov@weizmann.ac.il

<https://spindynamics.org/wiki/index.php?title=sim2liouv.m>
