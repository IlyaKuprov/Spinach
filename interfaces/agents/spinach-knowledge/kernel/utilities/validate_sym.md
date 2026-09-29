# kernel/utilities/validate_sym.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/validate_sym.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/validate_sym.m)

## Purpose

Extended validation of user-declared permutation symmetry. Confirms that the Zeeman, coupling, and giant-spin interactions stored in the spin system object are strictly invariant under every operation of each declared permutation group, so that the irreducible representation projectors built by `symmetry.m` correspond to a symmetry that the interactions actually possess.

The declared symmetry is a permutation of spin labels, not a spatial rotation; interaction tensors related by a rotation rather than being identical are not accepted. When no symmetry is declared, the function returns without performing any checks.

## Behaviour

- Syntax: `validate_sym(spin_system,bas)`.
- Calls `grumble(spin_system,bas)` first to enforce consistency of the basis specification fields.
- Returns immediately if `'symmetry'` is listed in `spin_system.sys.disable`.
- Returns immediately if `bas.sym_group` is absent or empty.
- Interaction agreement tolerance is `tol = 2*pi*spin_system.tols.inter_cutoff`, in rad/s.
- For each declared symmetry group `m`:
  - Retrieves the spin index vector `spins = bas.sym_spins{m}`.
  - Errors if the spins in the group are not all the same isotope.
  - Obtains the permutation elements of the declared group via `perm_group(bas.sym_group{m})`.
  - For each group operation `n = 1:group.order`, builds the global spin permutation `perm = 1:nspins` with `perm(spins) = spins(group.elements(n,:))`.
  - For each spin `k` in the symmetry group:
    - Checks Zeeman tensor invariance: errors if `norm(zeeman{perm(spins(k))}-zeeman{spins(k)},2) > tol`.
    - Checks giant-spin coefficient invariance across all ranks: errors if the coefficient cell arrays differ in length or any rank's coefficients differ by more than `tol` in the 2-norm.
    - Checks coupling tensor invariance against every spin `q = 1:nspins`: errors if `norm(get_coupling(spin_system,perm(spins(k)),perm(q)) - get_coupling(spin_system,spins(k),q),2) > tol`. The condition enforced is componentwise identity in the laboratory frame, because `symmetry.m` permutes spin labels without rotating them.
- On success, reports via `report(spin_system,'declared permutation symmetry is consistent with the interaction data.')`.
- The `grumble` subfunction enforces, when `bas.sym_group` is present:
  - `bas.sym_group` must be a cell array of strings.
  - `bas.sym_spins` must be specified alongside `bas.sym_group`.
  - `bas.sym_spins` must be a cell array of vectors.
  - `bas.sym_group` and `bas.sym_spins` must have the same number of elements.
  - Each `bas.sym_spins{n}` must contain only labels within `1..spin_system.comp.nspins` and at least 2 elements.

This is a service function of the Spinach kernel that should not be called directly; it is called by `symmetry.m`.

## Inputs and outputs

**Inputs:**

- `spin_system` — Spinach spin system description object, with the interaction arrays already processed by `create.m`.
- `bas` — basis specification structure; the fields used here are `bas.sym_group` (cell array of group names) and `bas.sym_spins` (cell array of spin index vectors), as described in the basis specification section of the online manual.

**Outputs:**

- None. The function throws a descriptive error when the interaction data does not obey a declared symmetry.

## References

- [validate_sym.m — Spinach Wiki](https://spindynamics.org/wiki/index.php?title=validate_sym.m)
