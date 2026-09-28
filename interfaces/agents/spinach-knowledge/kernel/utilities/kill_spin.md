# kernel/utilities/kill_spin.m

- Signature: `spin_system=kill_spin(spin_system,hit_list)`

## Purpose

Remove selected particles from a Spinach `spin_system` and update particle-indexed composition, interaction, relaxation, mode, and chemistry data.

## Physical / mathematical content

The routine edits the unified particle list and the associated data structures. It does not perform a spin-dynamics calculation.

## Numerical / algorithmic content

Selected entries are removed from particle-indexed arrays and matrices, and the remaining indices are compacted. Particle-indexed bosonic-mode parameters and pair channels are updated; nested modulation data are reindexed and entries involving removed particles are discarded. If no C, V, or T particle remains, the mode container is removed. Basis, connectivity, symmetry, and assumption data are invalidated or removed and must be rebuilt before subsequent calculations. Removing particles is rejected if it would leave fewer than two radical-pair electrons while radical-pair rates are present.

## Parameters / inputs

- `spin_system` - Spinach data structure whose particles are to be removed.
- `hit_list` - positive integer particle indices, or a logical mask with one element per particle.

## Outputs

- `spin_system` - updated data structure. Basis, connectivity, symmetry, and assumption information is not retained as valid data.

## Implementation structure

The function first validates the index list or mask, then removes the selected entries from composition and particle-indexed fields. It updates mode, relaxation, kinetics, and radical-pair bookkeeping, and removes derived basis/connectivity/symmetry/assumption information that is no longer valid.
