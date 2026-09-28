# kernel/utilities/shift_iso.m

- Signature: `tensors=shift_iso(tensors,spin_numbers,new_iso)`

## Purpose

Replaces the isotropic parts of selected interaction tensors with user-supplied values while retaining their anisotropic parts. This is useful for correcting DFT calculations when the anisotropy of the spin interactions is satisfactory but the isotropic part is not.

## Parameters / inputs

- `tensors` — cell array of interaction tensors as real 3×3 matrices.
- `spin_numbers` — indices of the tensors whose isotropic parts should be replaced; each index must be a positive integer within the cell array.
- `new_iso` — real replacement isotropic values, in the same order as `spin_numbers`. The two inputs must have the same number of elements.

## Outputs

- `tensors` — cell array of interaction tensors as 3×3 matrices, with the selected isotropic parts replaced.

## Implementation

For each index in `spin_numbers`, the function uses `mat2sphten` to separate the selected tensor into isotropic, rank-1, and rank-2 parts. It discards the original isotropic part and rebuilds the tensor as `sphten2mat([],rank1,rank2)+new_iso(n)*eye(3)`. Tensors not selected by `spin_numbers` are left unchanged. Inputs are checked for the stated types, dimensions, index bounds, and matching numbers of indices and replacement values.

[Source page](https://spindynamics.org/wiki/index.php?title=shift_iso.m)

Contacts: ledwards@cbs.mpg.de; ilya.kuprov@weizmann.ac.il.