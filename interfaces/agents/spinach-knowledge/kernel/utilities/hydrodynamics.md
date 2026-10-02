# kernel/utilities/hydrodynamics.m

## Purpose

A basic hydrodynamics infrastructure provider that returns first derivative operators with respect to the three sample coordinates ([source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/hydrodynamics.m)).

## Behaviour

- Syntax: `[Fx,Fy,Fz]=hydrodynamics(spin_system,parameters)`.
- Input consistency is enforced by an internal `grumble` subfunction, which validates `parameters.dims`, `parameters.npts`, and `parameters.deriv` and errors out with specific messages on invalid input.
- Derivative operators are built according to `parameters.deriv{1}`:
  - `'period'`: n-point central finite-difference matrices with periodic boundary conditions, obtained from `fdmat(...)` and divided by the grid spacing `parameters.dims(k)/parameters.npts(k)` for each dimension.
  - `'fourier'`: Fourier differentiation matrices from `fourdif(spin_system,...)`, scaled by `2*pi/parameters.dims(k)` for each dimension; with polyadics enabled, `fourdif` supplies the same action without storing dense derivative cores.
  - Any other value raises the error `'unrecognized derivative operator type.'`.
- The 1D derivative matrices are combined into full-space operators via `polyadic`, with identity matrices from `opium(...)` in the other dimensions, and multiplied by `-1i`:
  - 1D: `Fx=-1i*polyadic({{Dx}})`.
  - 2D: `Fx=-1i*polyadic({{opium(parameters.npts(2),1),Dx}})` and `Fy=-1i*polyadic({{Dy,opium(parameters.npts(1),1)}})`.
  - 3D: `Fx=-1i*polyadic({{opium(parameters.npts(3),1),opium(parameters.npts(2),1),Dx}})`, `Fy=-1i*polyadic({{opium(parameters.npts(3),1),Dy,opium(parameters.npts(1),1)}})`, and `Fz=-1i*polyadic({{Dz,opium(parameters.npts(2),1),opium(parameters.npts(1),1)}})`.
- The direct product order is Z(x)Y(x)X(x)Spin, corresponding to a column-wise vectorisation of a 3D array with dimensions ordered as [X Y Z].
- With Fourier derivatives and polyadics enabled, the outputs are action-only: neither `inflate()` nor `full()` can materialise their implicit FFT cores. Disable polyadics to obtain sparse numeric operators. Finite-difference polyadics remain materialisable.

## Inputs and outputs

Inputs:

- `spin_system` — spin system structure (used only to check `spin_system.sys.enable` for the `'polyadic'` option).
- `parameters.dims` — dimensions of the sample (meters), a one-, two-, or three-element row vector of positive real numbers.
- `parameters.npts` — number of points in each dimension of the sample, a one-, two-, or three-element row vector of real integers; must have the same number of elements as `parameters.dims`, and every entry must be at least 10 (fewer than 10 points triggers the error `'a spatial grid with fewer than 10 points is not a good idea - use a bigger grid.'`).
- `parameters.deriv` — a cell array with one or two elements: `{'fourier'}` requests Fourier differentiation matrices; `{'period',n}` requests n-point central finite-difference matrices with periodic boundary conditions. For `'period'`, the stencil size `n` must be a positive integer no greater than 7 (a larger value triggers the error `'differentiation stencil size greater than 7 is not a good idea - use a bigger grid.'`); for `'fourier'`, no second element is allowed.

Outputs:

- `Fx`, `Fy`, `Fz` — derivative operators, SI units. Outputs for dimensions absent from the grid remain empty (`[]`).

## References

- Source code: [kernel/utilities/hydrodynamics.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/hydrodynamics.m)
- Spinach Wiki: [hydrodynamics.m](https://spindynamics.org/wiki/index.php?title=hydrodynamics.m)

The implicit Fourier route preserves the first-derivative zero Nyquist convention on even grids. In `v2fplanck` and `imaging`, flow and diffusion inherit these actions; diffusion remains composed from the original first-derivative products rather than silently substituting a different second derivative. The finite-difference route is unchanged. Upload a resulting polyadic once with `gpuArray` when GPU actions are needed.
