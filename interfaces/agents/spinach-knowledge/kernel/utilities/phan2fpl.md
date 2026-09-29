# kernel/utilities/phan2fpl.m

## Purpose

Projects a spatial intensity distribution (a phantom) into Fokker–Planck space, using the phantom as the image painted by the supplied spin state. The function is documented as part of the Spinach library and referenced from the Spin Dynamics Wiki.

## Behaviour

- The function first validates its inputs via an internal consistency check (`grumble`).
- The phantom array is stretched into a column vector with `phan(:)` and combined with the spin state vector using the Kronecker product: `rho = kron(phan(:), rho)`.
- The result is a Fokker–Planck space state vector in which each spatial voxel of the phantom is paired with the corresponding spin state.
- Validation rules enforced by `grumble`:
  - `rho` must be numeric and a column vector (`size(rho,2) == 1`), otherwise an error is raised.
  - `phan` must be numeric, real, and 1D, 2D, or 3D (`ndims` in `[1 2 3]`), otherwise an error is raised.

## Inputs and outputs

Syntax:

```matlab
rho = phan2fpl(phan, rho)
```

Inputs:

- `phan` — phantom; the spatial distribution of the amplitude of the specified spin state. Must be a real numeric 1D, 2D, or 3D array.
- `rho` — Liouville space state vector. Must be a numeric column vector.

Outputs:

- `rho` — Fokker–Planck state vector.

## References

- Spin Dynamics Wiki page for `phan2fpl.m`: <https://spindynamics.org/wiki/index.php?title=phan2fpl.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/phan2fpl.m>
