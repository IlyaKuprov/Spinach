# kernel/utilities/kill_spin.m

- Signature: `spin_system=kill_spin(spin_system,hit_list)`

## Purpose

Removes the specified particles (spins or bosonic modes) from the spin_system structure and updates it accordingly. Syntax: spin_system=kill_spin(spin_system,hit_list)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

Retained bosonic-mode frequencies, carriers, anharmonicities, damping, and dephasing follow the new particle numbering. Pair channels and mode-pair modulation entries are reindexed on both particle axes; the nested spin tensor and field derivative blocks are also reindexed, with interactions involving removed particles discarded. Mode assumption strengths are cleared even when modes remain. When no C, V, or T particles remain, the mode container is discarded, permitting ordinary spin-only assumptions and retention options. Rebuild the basis and assumptions after removal.

## Parameters / inputs

- spin_system -primary Spinach data structure
- hit_list -a vector of integers or a logical
- vector giving the particle numbers
- to be removed from the system

## Outputs

- spin_system - the data structure with the indicated particles (spins or bosonic modes) and dependent information (basis, assumptions) removed
- Notes: basis, connectivity, symmetry, and assumption information
- is destroyed by this function; you would need to call the
- basis.m and assume.m functions again.

## Header notes

hit_list selects unified particle indices by an integer or logical vector. Removing particles invalidates basis and assumption information; rebuild both before subsequent calculations.
