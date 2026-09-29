# etc/molecules/fatty_acid.m

- MATLAB implementation: [etc/molecules/fatty_acid.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/molecules/fatty_acid.m)

**Call:** `[sys, inter] = fatty_acid(nprotons)` (for example, `[sys,inter]=fatty_acid(9)`).

**Input:** `nprotons` must be a finite, real, scalar, odd integer of at least 5. It sets the total number of spins; the routine rejects even values, nonintegers, nonfinite values, and values below 5.

## Model and numerical choices

This is a schematic all-proton fatty-acid-chain spin system, not a geometry-derived or compound-specific parameter set. The first three protons represent a methyl group and all receive chemical shift 1.0. The remaining proton shifts are assigned in alternating index groups (4,6,… and 5,7,…), each following a linear ramp from 1.5 to 5.0 across `(nprotons-3)/2` positions.

The first three spins form a methyl coupling triangle with each pair assigned 20. Their couplings to the first two chain positions (spins 4 and 5) are each assigned 7. Along the chain, the code generates selected scalar couplings using `15 + randn()` for one alternating set of adjacent pairs; the other adjacent set and selected next-neighbour/longer pairs use `7 + randn()`. Because `randn` is called during construction and the function does not set a random seed, those coupling entries can differ between calls; the output is an approximate example, not a reproducible fixed parameter table.

## Output

- `sys`: `nprotons` entries, all isotope `'1H'`.
- `inter`: chemical-shift and scalar-coupling cells populated by the assignments above. The function does not provide coordinates or identify a particular fatty-acid molecule.

**Source reference:** [Spinach Wiki: fatty_acid.m](https://spindynamics.org/wiki/index.php?title=fatty_acid.m).
