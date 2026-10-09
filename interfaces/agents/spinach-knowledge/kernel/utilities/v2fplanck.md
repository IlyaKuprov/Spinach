# kernel/utilities/v2fplanck.m

## Purpose

`v2fplanck.m` translates a stationary 3D velocity field and a diffusion tensor field into a Fokker-Planck evolution generator for use in Spinach spin dynamics simulations.

## Mathematical action

The spatial generator represents advection and diffusion on one to three voxel axes. A uniform velocity component multiplies its translation generator `F_i`; a spatially varying component contributes `diag(F_i*u_i)+diag(u_i)*F_i`, including the velocity-gradient term. Missing velocity components are zero. An isotropic diffusion coefficient adds `-1i*D*F_i*F_i` for each active direction; a constant anisotropic tensor adds cross-direction terms, while voxel-dependent diffusion uses `-1i*F_i*diag(D_ij)*F_j`. These are spatial transport terms, not spin Hamiltonian terms.

The spatial operator is tensored with the spin identity of dimension `bas.offsets(end)`, covering the full substance direct sum; in three dimensions its direct-product ordering is Z⊗Y⊗X⊗Spin. `F` is a sparse numeric matrix by default. Only when `spin_system.sys.enable` contains `'polyadic'` are the spatial derivative operators kept polyadic and the result polyadic. With Fourier derivatives, that result has implicit FFT cores and is action-only: `inflate()` cannot materialise it. Disable polyadics for a sparse numeric generator. Finite-difference polyadics remain materialisable.

## Valid parameter domain

`parameters.npts` contains one to three integer voxel counts, each at least 10; `parameters.dims` supplies the same number of positive spatial extents. The derivative scheme is `fourier` or periodic finite differences with a positive stencil size no greater than 7. Constant `parameters.diff` (nonnegative scalar or symmetric positive-semidefinite tensor) and voxel-wise tensor components are alternative specifications, not simultaneous inputs; a 3D voxel-wise tensor requires all nine Cartesian components.

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

- `F` — the spatial Fokker–Planck generator, sparse numeric by default or polyadic when that feature is enabled.

## References

- Source code: [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/v2fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/v2fplanck.m)
- Spinach Wiki: [https://spindynamics.org/wiki/index.php?title=v2fplanck.m](https://spindynamics.org/wiki/index.php?title=v2fplanck.m)
