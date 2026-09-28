# kernel/utilities/v2fplanck.m

- Signature: `F=v2fplanck(spin_system,parameters)`

## Purpose

Translates a stationary 3D velocity field and a diffusion tensor field into a Fokker–Planck evolution generator.

## Physical / mathematical content

The generator combines spatial transport and diffusion with the spin system. Spatial coordinates are represented on a 3D grid.

## Numerical / algorithmic content

The function obtains the spatial generators from `hydrodynamics` and returns the result as a polyadic object. Use `inflate()` to obtain the corresponding sparse matrix.

## Parameters / inputs

- `parameters.u` — X components of the velocity vectors for each voxel in the sample, m/s; a scalar specifies spatially uniform flow along X.
- `parameters.v` — Y components of the velocity vectors for each voxel in the sample, m/s; a scalar specifies spatially uniform flow along Y.
- `parameters.w` — Z components of the velocity vectors for each voxel in the sample, m/s; a scalar specifies spatially uniform flow along Z.
- `parameters.diff` — diffusion coefficient or 3×3 tensor, m^2/s, for situations where this parameter is the same in every voxel.
- `parameters.dxx`, `parameters.dxy`, …, `parameters.dzz` — Cartesian components of the diffusion tensor for each voxel of the sample.
- `parameters.dims` — dimensions of the 3D box, meters.
- `parameters.npts` — number of points in each dimension of the 3D box.
- `parameters.deriv` — `{'fourier'}` uses Fourier differentiation matrices; `{'period',n}` requests n-point central finite-difference matrices with periodic boundary conditions.

## Outputs

- `F` — spatial dynamics generator.

- Note: the direct product order is `Z(x)Y(x)X(x)Spin`; this corresponds to a column-wise vectorization of a 3D array with dimensions ordered as `[X Y Z]`.
- Note: polyadic objects are returned; use `inflate()` to get the corresponding sparse matrix.

## Implementation structure

- Checks the parameter structure with `grumble` and obtains spatial generators from `hydrodynamics`.
- Constructs the Fokker–Planck generator using the velocity and diffusion inputs.
- Source documentation: <https://spindynamics.org/wiki/index.php?title=v2fplanck.m>
