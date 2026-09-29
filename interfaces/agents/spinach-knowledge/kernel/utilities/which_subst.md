# kernel/utilities/which_subst.m

## Purpose

Determines which substance in a spin system hosts a specified list of spins, throwing an error if the spins span more than one substance or belong to none. Source: [Spinach GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/which_subst.m).

## Behaviour

- Syntax: `subst=which_subst(spin_system,spins)`.
- Validates the spin list via an internal `grumble` subfunction: spins must be a real numeric vector of positive integers, must not exceed `spin_system.comp.nspins`, and must contain no repeated entries.
- Builds a logical mask over `spin_system.chem.parts`, marking substances whose spin list contains any of the specified spins.
- Errors with `'spin list crosses chemical boundaries.'` if more than one substance matches, or `'spins do not belong to any substance.'` if none matches.
- Returns the index of the single matching substance via `find` on the mask.
- Performs a final check that all specified spins are members of `spin_system.chem.parts{subst}`; otherwise errors again with `'spin list crosses chemical boundaries.'`.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system structure containing `chem.parts` (cell array of per-substance spin lists) and `comp.nspins` (total spin count).
- `spins` — a list of positive integers identifying spins.

**Outputs**

- `subst` — a positive integer giving the substance number hosting all specified spins.

## References

- [which_subst.m — Spinach Wiki](https://spindynamics.org/wiki/index.php?title=which_subst.m)
- [which_subst.m — GitHub source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/which_subst.m)
