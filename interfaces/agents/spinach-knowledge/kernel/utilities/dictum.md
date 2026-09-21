# kernel/utilities/dictum.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/dictum.m`
- Signature: `spin_system=dictum(spin_system,spins,strength)`
- Total lines: 153

## Purpose

Overrides default assumptions about interaction terms surviving rotating frame transformations. Syntax: spin_system=dictum(spin_system,spins,strength)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

- Overrides default assumptions about interaction terms surviving
- rotating frame transformations. Syntax:
- spin_system=dictum(spin_system,spins,strength)
- spin_system -Spinach spin system information
- object coming out of assume.m
- spins -a vector with one or two numbers
- or a cell array with one or two
- strings, e.g. [2 4] or {'1H'},
- where one element would cause
- Zeeman interaction assumptions
- to be modified, and two elements
- would cause coupling assumptions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `spins()`, `isscalar()`, `iscell()`, `cellfun()`, `strcmp()`, `isfield()`, `assume()`, `ischar()`, `isvector()`, `any()`, `ismember()`.
