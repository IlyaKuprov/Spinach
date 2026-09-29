# kernel/utilities/hydrodynamics.m

## Purpose

A basic hydrodynamics infrastructure provider that returns first derivative operators with respect to the three sample coordinates ([source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/hydrodynamics.m)).

## Behaviour

- Syntax: `[Fx,Fy,Fz]=hydrodynamics(spin_system,parameters)`.
- Input consistency is enforced by an internal `grumble` subfunction, which validates `parameters.dims`, `parameters.npts`, and `parameters.deriv` and errors out with specific messages on invalid input.
- Derivative operators are built according to `parameters.deriv{1}`:
  - `'period'`: n-point central finite-difference matrices with periodic boundary conditions, obtained from `fdmat(...)` and divided by the grid spacing `parameters.dims(k)/parameters.npts(k)` for each dimension.
  - `'fourier'`: Fourier differentiation matrices from `fourdif(...)`, scaled by `2*pi/parameters.dims(k)` for each dimension.
  - Any other value raises the error `'unrecognized derivative operator type.'`.
- The 1D derivative matrices are combined into full-space operators via `polyadic`, with identity matrices from `opium(...)` in the other dimensions, and multiplied by `-1i`:
  - 1D: `Fx=-1i*polyadic({{Dx}})`.
  - 2D: `Fx=-1i*polyadic({{opium(parameters.npts(2),1),Dx}})` and `Fy=-1i*polyadic({{Dy,opium(parameters.npts(1),1)}})`.
  - 3D: `Fx=-1i*polyadic({{opium(parameters.npts(3),1),opium(parameters.npts(2),1),Dx}})`, `Fy=-1i*polyadic({{opium(parameters.npts(3),1),Dy,opium(parameters.npts(1),1)}})`, and `Fz=-1i*polyadic({{Dz,opium(parameters.npts(2),1),opium(parameters.npts(1),1)}})`.
- The direct product order is Z(x)Y(x)X(x)Spin, corresponding to a column-wise vectorisation of a 3D array with dimensions ordered as [X Y Z].
- Polyadic objects are returned; unless `'polyadic'` is listed in `spin_system.sys.enable`, the outputs are passed through `inflate()` to yield the corresponding sparse matrices.

## Inputs and outputs

Inputs:

- `spin_system` — spin system structure (used only to check `spin_system.sys.enable` for the `'polyadic'` option).
- `parameters.dims` — dimensions of the sample (meters), a one-, two-, or three-element row vector of positive real numbers.
- `parameters.npts` — number of points in each dimension of the sample, a one-, two-, or three-element row vector of real integers; must have the same number of elements as `parameters.dims`, and every entry must be at least 10 (fewer than 10 points triggers the error `'a spatial grid with fewer than 10 points is not a good idea - use a bigger grid.'`).
- `parameters.deriv` — a cell array with one or two elements: `{'fourier'}` requests Fourier differentiation matrices; `{'period',n}` requests n-point central finite-difference matrices with periodic boundary conditions. For `'period'`, the stencil size `n` must be a positive integer no greater than 7 (a larger value triggers the error `'differentiation stencil size greater than 7 is not a good idea - use a bigger grid.'`); for `'fourier'`, no second element is allowed.

Outputs:

- `Fx`, `Fy`, `Fz` — derivative matrices, SI units. Outputs for dimensions absent from the grid remain empty (`[]`).

## References

- Source code: [kernel/utilities/hydrodynamics.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/hydrodynamics.m)
- Spinach Wiki: [hydrodynamics.m](https://spindynamics.org/wiki/index.php?title=hydrodynamics.m)
