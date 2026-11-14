# kernel/cce/cce_bath.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_bath.m)

## Purpose

Construct explicit, lattice or random spin baths with explicit isotope assignments.

## Behaviour

- Builds Spinach sys/inter inputs and a physical bath record; cell vectors and coordinates use Angstrom, magnetic field uses tesla.
- Supports seeded random coordinates/occupancies and exact or sampled occupancy configurations. Missing occupancy probability denotes vacancy.
- Does not add point dipoles. Optional tensors_Hz and pairs_Hz carry Cartesian tensor data; rotation, translation_A, remove_ids, selection and groups retain stable site IDs.

## Inputs and outputs

```matlab
[sys,inter,bath]=cce_bath(specification)
```

**Inputs**

- `specification` — mode, isotope cell row, magnetic field and mode-specific coordinates, cell/fractional/repeats, or bounds_A/count-or-density/dimension/seed; optional occupancy and interaction declarations.

**Outputs**

- `sys`, `inter` — inputs for Spinach create.
- `bath` — coordinates, physical identity, occupancy/configuration and construction metadata.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_bath.m)
- [cce_import.m](cce_import.md)
- [cce_prepare.m](cce_prepare.md)
