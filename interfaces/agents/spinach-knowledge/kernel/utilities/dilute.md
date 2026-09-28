# kernel/utilities/dilute.m

- Signature: `subsystems=dilute(spin_system,isotope,tuples)`

## Purpose

Splits a spin system into subsystems, each retaining a selected tuple of spins of the specified dilute isotope. Other spins of that isotope are removed from each subsystem, and spin system data is updated accordingly. Basis set information, if present, is destroyed.

## Physical / mathematical content

- Each returned subsystem represents one selection of `tuples` instances of the specified isotope.

## Numerical / algorithmic content

- If `tuples` is omitted, it defaults to `1`.
- The function finds spins whose isotope string matches `isotope`, then generates the distinct combinations of `tuples` such spins.
- For each combination, `kill_spin` removes the other spins of that isotope from a copy of the spin system.
- `isotope` must be a character string, `tuples` must be a positive integer, and the system must contain at least `tuples` instances of the isotope.

## Syntax

```matlab
subsystems=dilute(spin_system,isotope,tuples)
```

## Parameters / inputs

- spin_system -spin system object obtained
- from create.m function
- isotope -isotope specification string,
- for example '13C'
- tuples -1 returns all subsystems with
- a single instance of the isoto-
- pe, 2 returns all subsystems
- where two instances of the iso-
- tope are present, etc.

## Outputs

- subsystems -a cell array of spin system
- objects, with each cell cor-
- responding to one of the new-
- ly formed isotopomers.

## Implementation structure

- Returns a cell array with one spin system per distinct combination of dilute-isotope spins.