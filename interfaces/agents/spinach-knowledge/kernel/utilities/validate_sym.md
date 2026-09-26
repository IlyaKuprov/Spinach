# kernel/utilities/validate_sym.m

- Signature: `validate_sym(spin_system,bas)`

## Purpose

Extended validation of user-declared permutation symmetry. Confirms that the Zeeman, coupling, and giant-spin interactions stored in the spin system object are strictly invariant under every operation of each declared permutation group, so that the irreducible representa- tion projectors built by symmetry.m correspond to a symmetry that the interactions actually possess. The declared symmetry is a permu- tation of sp

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach spin system description object, with the
- interaction arrays already processed by create.m
- bas -basis specification structure; the fields used
- here are bas.sym_group (cell array of group names)
- and bas.sym_spins (cell array of spin index vec-
- tors), as described in the basis specification
- section of the online manual

## Outputs

- this function returns nothing; it throws a descriptive error when
- the interaction data does not obey a declared symmetry
- Note: this is a service function of the Spinach kernel that should
- not be called directly; it is called by symmetry.m

## Implementation structure

- Extended validation of user-declared permutation symmetry. Confirms
- that the Zeeman, coupling, and giant-spin interactions stored in the
- spin system object are strictly invariant under every operation of
- each declared permutation group, so that the irreducible representa-
- tion projectors built by symmetry.m correspond to a symmetry that
- the interactions actually possess. The declared symmetry is a permu-
- tation of spin labels, not a spatial rotation; interaction tensors
- related by a rotation rather than being identical are not accepted.
- When no symmetry is declared, the function returns without perfor-
- ming any checks. Syntax:
- validate_sym(spin_system,bas)
- spin_system -Spinach spin system description object, with the
