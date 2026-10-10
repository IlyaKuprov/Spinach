# kernel/utilities/shift_iso.m

## Purpose

Replaces the isotropic parts of interaction tensors with user-supplied values. This is useful for correcting DFT calculations, where the anisotropy of the various spin interactions is usually satisfactory, but the isotropic part is not.

Source: [kernel/utilities/shift_iso.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/shift_iso.m)

## Behaviour

- Syntax: `tensors=shift_iso(tensors,spin_numbers,new_iso)`.
- The function first runs a consistency check (`grumble`) on all inputs.
- For each spin number in `spin_numbers`, the corresponding tensor is decomposed into spherical tensor components using `mat2sphten`, which returns the isotropic part (`rank0`, discarded), the rank-1 component, and the rank-2 component.
- The tensor is then rebuilt with `sphten2mat([],rank1,rank2)`, keeping only the anisotropic parts, and the new isotropic contribution is added as `new_iso(n)*eye(3)`.
- Only the tensors listed in `spin_numbers` are modified; all other tensors in the cell array are returned unchanged.

## Inputs and outputs

**Inputs**

- `tensors` — a cell array of interaction tensors as 3x3 matrices.
- `spin_numbers` — a vector containing the numbers of spins in the `tensors` array that should have the isotropic parts replaced. Must be a vector of positive integers, each no greater than the number of tensors supplied.
- `new_iso` — a vector containing the new isotropic parts in the same order as the spin numbers listed in `spin_numbers`. Must be a vector of real numbers with the same number of elements as `spin_numbers`.

**Outputs**

- `tensors` — a cell array of interaction tensors as 3x3 matrices, with the isotropic parts of the selected tensors replaced.

**Validation errors**

- `tensors` must be a cell array of real 3x3 matrices (empty elements are permitted).
- `spin_numbers` must be a vector of positive integers.
- No index in `spin_numbers` may exceed the number of tensors supplied.
- `new_iso` must be a vector of real numbers.
- `spin_numbers` and `new_iso` must have the same number of elements.

## References

- Spinach documentation: [shift_iso.m](https://spindynamics.org/wiki/index.php?title=shift_iso.m)
- Source code: [kernel/utilities/shift_iso.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/shift_iso.m)
