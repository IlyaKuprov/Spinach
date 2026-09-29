# kernel/utilities/dictum.m

## Purpose

Overrides default assumptions about which interaction terms survive rotating frame transformations, as originally set by `assume.m`. One spin selection modifies Zeeman interaction assumptions; a two-spin selection modifies coupling assumptions.

## Behaviour

- Syntax: `spin_system=dictum(spin_system,spins,strength)`.
- Calls the internal `grumble` function first to enforce consistency: `spin_system.inter.coupling.strength` and `spin_system.inter.zeeman.strength` fields must exist (i.e. `assume()` must have been run before calling this function), `strength` must be a character string, and `spins` must be either a vector of one or two positive integers not exceeding `spin_system.comp.nspins`, or a cell array of one or two character strings matching isotopes present in the system. Otherwise the function errors out.
- Numerical specification with two spins: reports the previous coupling assumption for the spin pair, sets `spin_system.inter.coupling.strength{spins(1),spins(2)}` and the transposed element to the new strength, and reports the new assumption.
- Numerical specification with one spin: reports the previous Zeeman assumption, sets `spin_system.inter.zeeman.strength{spins}` to the new strength, and reports the new assumption.
- Isotope specification with two strings: loops over all spin pairs whose isotope names match the two given strings, and for each pair reports the previous coupling assumption, sets both symmetric `spin_system.inter.coupling.strength{n,k}` and `{k,n}` entries, and reports the new assumption.
- Isotope specification with one string: loops over all spins whose isotope name matches, and for each reports the previous Zeeman assumption, sets `spin_system.inter.zeeman.strength{n}`, and reports the new assumption.
- Any other input syntax raises the error `'incorrect input syntax.'`.
- Reporting is done via the `report` function, printing the spin indices, isotope names, the old assumption, and the new assumption.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach spin system information object coming out of `assume.m`.
- `spins` — a vector with one or two numbers, or a cell array with one or two strings, e.g. `[2 4]` or `{'1H'}`. One element causes Zeeman interaction assumptions to be modified; two elements cause coupling assumptions to be modified.
- `strength` — new strength specification; see the source code of `assume.m` for the available strength specs.

**Outputs**

- `spin_system` — updated Spinach spin system information object that will be used by `hamiltonian.m` to build the Hamiltonian.

## References

- Source: [kernel/utilities/dictum.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dictum.m)
- Spin Dynamics Wiki: [dictum.m](https://spindynamics.org/wiki/index.php?title=dictum.m)
