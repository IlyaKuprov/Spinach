# kernel/utilities/sim2liouv.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sim2liouv.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sim2liouv.m)

## Purpose

Moves a zeeman-hilb simulation context into Liouville space. When the formalism specified in the spin system object is `zeeman-hilb`, this function projects the evolution generators into Liouville space, converts the standard state-like and operator-like fields of the parameters structure, rebuilds the basis index table, migrates the symmetry irrep projectors into the adjoint representation, and sets the formalism to `zeeman-liouv`. For all other formalisms, every argument is returned unchanged. This makes Liouville-space pulse sequences callable with zeeman-hilb inputs.

## Behaviour

- Calls `grumble` to enforce consistency: the formalism in `spin_system.bas.formalism` must be one of `sphten-liouv`, `zeeman-liouv`, `zeeman-hilb`, or `zeeman-wavef`; `parameters` must be a structure; `H`, `R`, and `K` must be numeric arrays.
- Only proceeds when the formalism is `zeeman-hilb`; otherwise all arguments are returned unchanged.
- Projects the evolution generators: `H` via `hilb2liouv(H,'comm')`, and `R` and `K` via `hilb2liouv(...,'acomm')`. Empty matrices are passed through.
- Stretches the state-like parameter fields `rho0`, `coil`, and `screen` (matrices or their horizontal concatenations) into state vectors using `reshape(...,hdim^2,[])`, where `hdim` is the Hilbert space dimension.
- Converts the operator-like parameter fields `pulse_op`, `mw_oper`, `ez_oper`, and `homodec_oper`, when present, into commutation superoperators via `hilb2liouv(...,'comm')`.
- Rebuilds the basis index table for Liouville space as `[repmat(zbas,[hdim 1]) kron(zbas,ones(hdim,1))]`, where `zbas` is the Hilbert space basis table.
- If `spin_system.bas.basis_hash` exists, refreshes it with `md5_hash` of the new basis.
- If `spin_system.bas.irrep` exists, migrates the Hilbert space irreps into the adjoint representation: each ordered pair of Hilbert space irrep projectors yields the Liouville space irrep projector `kron(conj(S(k)),S(n))` with dimension the product of the pair dimensions. Diagonal irrep pairs, which share the unit state, are merged into the first subspace; off-diagonal pairs are stored as separate subspaces. The resulting number of subspaces is `n_irreps^2-n_irreps+1`.
- Sets `spin_system.bas.formalism` to `zeeman-liouv`.
- If `R` is non-empty, projects the unit state out of the relaxation superoperator using `R=R-U*(U'*R)-(R*U)*U'+U*(U'*R*U)*U'`, where `U=unit_state(spin_system)`, so the unit state is neither damped nor a source of relaxation and the trace is conserved.
- A Hilbert space density matrix block `S(n)*Y*S(k)'` maps to the state vector `kron(conj(S(k)),S(n))*Y(:)`. Every such subspace is invariant under superoperators built from symmetry-respecting Hilbert space generators; unpopulated subspaces are dropped by `reduce.m` at run time.
- The anticommutation superoperator of a Hilbert space relaxation matrix damps the unit state, which the Liouville space branch of `relaxation.m` never does. For the scalar damping matrix built in Hilbert space by the kernel, the projected result coincides with the Liouville space damp operator. The unit state spans all diagonal irrep pairs, which are merged into one subspace so the projected `R` stays block-diagonal in the irrep table that `reduce.m` evolves independently.
- Reports progress messages via `report` when projecting, migrating irreps, and exempting the unit state.

## Inputs and outputs

**Syntax:**

```
[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)
```

**Inputs:**

- `spin_system` — Spinach spin system object.
- `parameters` — pulse sequence parameters structure; the state-like fields `rho0`, `coil`, and `screen` (matrices or their horizontal concatenations) are stretched into state vectors, and the operator-like fields `pulse_op`, `mw_oper`, `ez_oper`, and `homodec_oper` become commutation superoperators, when present.
- `H` — Hamiltonian operator, converted into a commutation superoperator; an empty matrix is passed through.
- `R` — relaxation matrix, converted into an anticommutation superoperator with the unit state exempted from damping; an empty matrix is passed through.
- `K` — kinetics matrix, converted into an anticommutation superoperator; an empty matrix is passed through.

**Outputs:**

- `spin_system` — spin system object with `zeeman-liouv` formalism and basis information.
- `parameters` — parameters structure with the standard fields converted into Liouville space.
- `H`, `R`, `K` — Liouville space evolution generators.

## References

- Spinach Wiki: [sim2liouv.m](https://spindynamics.org/wiki/index.php?title=sim2liouv.m)
