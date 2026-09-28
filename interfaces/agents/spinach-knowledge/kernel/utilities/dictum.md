# kernel/utilities/dictum.m

- Signature: `spin_system=dictum(spin_system,spins,strength)`

## Purpose

Overrides default assumptions about which interaction terms survive rotating-frame transformations.

## Physical / mathematical content

One spin selects a Zeeman interaction assumption; two spins select a coupling assumption.

## Numerical / algorithmic content

Numeric spin indices update the specified spin or pair. Isotope strings select matching spins in the system: one isotope updates each matching Zeeman assumption, while two isotopes update coupling assumptions for matching pairs. Coupling updates are written in both index orders. The function reports the previous and new assumptions.

## Parameters / inputs

- spin_system -Spinach spin system information
- object coming out of assume.m
- spins -a vector with one or two numbers
- or a cell array with one or two
- strings, e.g. [2 4] or {'1H'},
- where one element would cause
- Zeeman interaction assumptions
- to be modified, and two elements
- would cause coupling assumptions
- to be modified.
- strength -new strength specification, see
- the source code of assume.m for
- the available strength specs

## Outputs

- spin_system -updated Spinach spin system in-
- formation object that will be
- used by hamiltonian.m to build
- the Hamiltonian

## Implementation structure

The function checks that `assume()` has supplied Zeeman and coupling strength information, that `strength` is a character string, and that `spins` contains one or two valid spin indices or isotope strings. It then updates the selected entries in `spin_system.inter.zeeman.strength` or `spin_system.inter.coupling.strength`.

Source reference: <https://spindynamics.org/wiki/index.php?title=dictum.m>