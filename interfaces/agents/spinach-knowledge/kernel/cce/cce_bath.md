# kernel/cce/cce_bath.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_bath.m)

## Purpose

Construct explicit, lattice or random spin baths with explicit isotope assignments.

## Behaviour

- `tensors_Hz` specifies total local Zeeman frequency tensors, not offsets: after rotation, nuclear tensors are converted to ppm and electron tensors to g units so `create` recovers `2*pi*tensors_Hz` in rad/s. A nonzero tensor at zero field is refused because native Zeeman inputs cannot represent it; zero tensors are accepted.

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

Joint isotope configurations contain one entry per pre-exclusion candidate, including vacancies. For random baths this length is the explicit count or the density-derived count, independently of whether the base isotope table contains a single label or one label per candidate. Site-resolved independent occupancy probabilities likewise use one row per random candidate.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_bath.m)
- [cce_import.m](cce_import.md)
- [cce_prepare.m](cce_prepare.md)

Bath constructors reject unknown fields, mutually exclusive random count/density, mode-inapplicable geometry, and incompatible explicit versus independent isotopologue specifications.

For random baths an explicit dimension, when supplied, must equal the number of nonzero Cartesian box extents. Selection subfields are center_A, inner_A and outer_A.
