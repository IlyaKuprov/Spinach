# kernel/utilities/kill_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/kill_spin.m`
- Signature: `spin_system=kill_spin(spin_system,hit_list)`
- Total lines: 194

## Purpose

Removes the specified spins from the spin_system structure and updates it accordingly. Syntax: spin_system=kill_spin(spin_system,hit_list)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -primary Spinach data structure
- hit_list -a vector of integers or a logical
- vector giving the numbers of spins
- to be removed from the system

## Outputs

- spin_system -the data structure with the indica-
- ted spins and all dependent infor-
- mation (basis, assumptions) removed
- Notes: basis, connectivity, symmetry, and assumption information
- is destroyed by this function; you would need to call the
- basis.m and assume.m functions again.

## Implementation structure

- Removes the specified spins from the spin_system structure and
- updates it accordingly. Syntax:
- spin_system=kill_spin(spin_system,hit_list)
- spin_system -primary Spinach data structure
- hit_list -a vector of integers or a logical
- vector giving the numbers of spins
- to be removed from the system
- spin_system -the data structure with the indica-
- ted spins and all dependent infor-
- mation (basis, assumptions) removed
- is destroyed by this function; you would need to call the
- basis.m and assume.m functions again.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `islogical()`, `report()`, `num2str()`, `isfield()`, `md5_hash()`, `srsk_spins()`, `false()`, `subsystem_idx()`, `true()`, `reacting_spins()`, `rmfield()`, `any()`.
