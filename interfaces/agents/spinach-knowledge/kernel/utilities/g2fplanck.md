# kernel/utilities/g2fplanck.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/g2fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/g2fplanck.m)

## Purpose

Returns gradient operators within the Fokker-Planck formalism used in the imaging module of Spinach.

## Behaviour

- Syntax: `G=g2fplanck(spin_system,parameters)`.
- Consistency is enforced by an internal `grumble` subfunction, which errors when: the primary magnet field is zero; `parameters.dims` is missing, non-numeric, non-real, non-positive, or has fewer than one or more than three elements; `parameters.npts` is missing, non-numeric, non-real, less than one, non-integer, or has fewer than one or more than three elements; or `parameters.dims` and `parameters.npts` differ in length.
- The magnet Zeeman Hamiltonian is built via `hamiltonian(assume(spin_system,'labframe','zeeman'))` and divided by `spin_system.inter.magnet`, giving a per-tesla Hamiltonian `H`.
- Gradients are assumed linear and centred on the middle of the sample; each spatial axis uses `linspace(-0.5,0.5,npts)` scaled by the corresponding box dimension, placed on the diagonal with `spdiags`.
- One dimension: `Gx` is the scaled diagonal operator combined with `H` via `polyadic({{Gx,H}})`.
- Two dimensions: `Gx` and `Gy` are built as `polyadic({{opium(npts(2),1),Gx,H}})` and `polyadic({{Gy,opium(npts(1),1),H}})`; if `parameters.grad_angles` is present, the two operators are rotated with a 2D rotation matrix built from `cos`/`sin` of the angle.
- Three dimensions: `Gx`, `Gy`, `Gz` are built as `polyadic({{opium(npts(3),1),opium(npts(2),1),Gx,H}})`, `polyadic({{opium(npts(3),1),Gy,opium(npts(1),1),H}})`, and `polyadic({{Gz,opium(npts(2),1),opium(npts(1),1),H}})`; if `parameters.grad_angles` is present, the operators are rotated using a direction cosine matrix from `euler2dcm`.
- Any other number of spatial dimensions raises the error `'incorrect number of spatial dimensions.'`.
- The direct product order is Z(x)Y(x)X(x)Spin, corresponding to a column-wise vectorisation of a 3D array with dimensions ordered as [X Y Z].
- Polyadic objects are returned; use `inflate()` to obtain the corresponding sparse matrix.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach spin system object; its `inter.magnet` field must be non-zero.
- `parameters.dims` — vector with one, two, or three elements giving the dimensions of the box, in metres.
- `parameters.npts` — vector with one, two, or three elements giving the number of points in each dimension of the box; must match the length of `parameters.dims`.
- `parameters.grad_angles` — optional; when present, gradient operators are rotated by the supplied angles (2D rotation matrix in the two-dimensional case, `euler2dcm` direction cosine matrix in the three-dimensional case).

**Outputs**

- `G` — cell array with the three gradient operators ordered as `{Gx,Gy,Gz}`, normalised to 1 T/m, with empty matrices for non-existent dimensions.

## References

- Spinach Wiki: [g2fplanck.m](https://spindynamics.org/wiki/index.php?title=g2fplanck.m)
