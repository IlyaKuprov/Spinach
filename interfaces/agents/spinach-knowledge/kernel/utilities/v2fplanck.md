# kernel/utilities/v2fplanck.m

## Purpose

`v2fplanck.m` translates a stationary 3D velocity field and a diffusion tensor field into a Fokker-Planck evolution generator for use in Spinach spin dynamics simulations.

## Behaviour

- The function is called as `F=v2fplanck(spin_system,parameters)`.
- It first validates the input parameters via an internal consistency-checking routine (`grumble`), which enforces:
  - `parameters.npts` must be a 1, 2, or 3-element row vector of positive integers, each at least 10.
  - Either `parameters.diff` or the voxel-wise `parameters.d**` tensor components may be specified, but not both.
  - For 1D samples, `parameters.v` and `parameters.w` are not applicable; only `parameters.u` and `parameters.dxx` are used.
  - For 2D samples, `parameters.w` is not applicable; `parameters.dxx`, `dxy`, `dyx`, `dyy` must all be specified simultaneously, the tensor field must be symmetric (within `1e-10` relative tolerance), and positive semidefinite (eigenvalues not below `-20*eps` scaled).
  - For 3D samples, all nine tensor components (`dxx` through `dzz`) must be specified simultaneously, with the same symmetry and positive semidefiniteness checks.
  - `parameters.diff`, if a matrix, must be symmetric with non-negative eigenvalues; if scalar, it must be non-negative.
  - `parameters.dims` must be a row vector of finite positive real numbers matching the length of `parameters.npts`.
  - `parameters.deriv` must be a cell array with first element `'fourier'` or `'period'`; if `'period'`, a second element specifying a positive integer stencil size no greater than 7 is required.
- The function obtains translation generators `Fx`, `Fy`, `Fz` from `hydrodynamics(spin_system,parameters)`.
- The number of voxels is computed as `prod(parameters.npts)`.
- Flow terms:
  - If a velocity component (`u`, `v`, or `w`) is a scalar, uniform flow is added as a scalar multiple of the corresponding translation generator.
  - If it is a vector/array, non-uniform flow is added using symmetric products of the translation generator with diagonal matrices built from the velocity field.
  - Velocity components default to 0 if not specified.
- Diffusion terms:
  - Isotropic scalar `diff` adds `-1i*diff*F*F` terms along each active dimension.
  - Anisotropic tensor `diff` adds cross terms such as `-1i*diff(1,2)*Fx*Fy` and `-1i*diff(2,1)*Fy*Fx` for all active dimension pairs.
  - Voxel-wise tensor components add terms of the form `-1i*Fx*spdiags(dxy(:),...)*Fy` for all active pairs.
- The generator is cleaned using `clean_up(spin_system,F,spin_system.tols.liouv_zero)`.
- Finally, the spatial generator is Kronecker-multiplied with the spin identity operator: `F=kron(F,opium(spn_dim,1))`, where `spn_dim` is the size of the spin basis.
- The direct product order is Z(x)Y(x)X(x)Spin, corresponding to column-wise vectorisation of a 3D array with dimensions ordered as [X Y Z].
- Polyadic objects are returned; use `inflate()` to obtain the corresponding sparse matrix.

## Inputs and outputs

**Inputs:**

- `spin_system` — the Spinach spin system object.
- `parameters` — a structure containing:
  - `parameters.u` — X components of velocity vectors for each voxel, m/s; a scalar specifies spatially uniform flow along X.
  - `parameters.v` — Y components of velocity vectors for each voxel, m/s; a scalar specifies spatially uniform flow along Y.
  - `parameters.w` — Z components of velocity vectors for each voxel, m/s; a scalar specifies spatially uniform flow along Z.
  - `parameters.diff` — diffusion coefficient or 3x3 tensor, m^2/s, for spatially uniform diffusion.
  - `parameters.dxx` ... `parameters.dzz` — Cartesian components of the diffusion tensor for each voxel.
  - `parameters.dims` — dimensions of the 3D box, meters.
  - `parameters.npts` — number of points in each dimension of the 3D box.
  - `parameters.deriv` — `{'fourier'}` uses Fourier differentiation matrices; `{'period',n}` requests n-point central finite-difference matrices with periodic boundary conditions.

**Outputs:**

- `F` — the spatial dynamics generator (Fokker-Planck evolution generator).

## References

- Source code: [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/v2fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/v2fplanck.m)
- Spinach Wiki: [https://spindynamics.org/wiki/index.php?title=v2fplanck.m](https://spindynamics.org/wiki/index.php?title=v2fplanck.m)
