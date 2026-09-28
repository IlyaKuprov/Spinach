# kernel/utilities/which_subst.m

- Signature: `subst=which_subst(spin_system,spins)`

## Purpose

Returns the number of the substance containing all specified spins. Raises an error if the spins cross chemical boundaries or do not belong to any substance.

## Parameters

- `spin_system`: spin system whose `chem.parts` lists the spins in each substance and whose `comp.nspins` gives the total spin count.
- `spins`: a list of distinct positive integer spin numbers, each no greater than the total spin count.

## Output

- `subst`: a positive integer identifying the substance.

## Behavior

The function validates `spins`, finds the substance containing them, and checks that every specified spin belongs to that substance. It raises an error for invalid, repeated, or out-of-range spin numbers; spins absent from all substances; or spins crossing chemical boundaries.

Source: <https://spindynamics.org/wiki/index.php?title=which_subst.m>