# kernel/utilities/hydrodynamics.m

- Signature: `[Fx,Fy,Fz]=hydrodynamics(spin_system,parameters)`

## Purpose

A basic hydrodynamics infrastructure provider, returns first derivative operators with respect to the three sample coordi- nates. Syntax: [Fx,Fy,Fz]=hydrodynamics(parameters)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Parameters / inputs

- parameters.dims -dimensions of the sample (meters),
- one, two, or three-element row
- vector
- parameters.npts -number of points in each dimension
- of the sample, one, two, or three-
- element row vector
- parameters.deriv -{'fourier'} requests Fourier diffe-
- rentiation matrices; {'period',n}
- requests n-point central finite-
- difference matrices with periodic
- boundary conditions

## Outputs

- Fx, Fy, Fz -derivative matrices, SI units
- Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
- responds to a column-wise vectorization of a 3D array
- with dimensions ordered as [X Y Z].
- Note: polyadic objects are returned, use inflate() to get the
- corresponding sparse matrix.

## Implementation structure

- A basic hydrodynamics infrastructure provider, returns first
- derivative operators with respect to the three sample coordi-
- nates. Syntax:
- [Fx,Fy,Fz]=hydrodynamics(parameters)
- parameters.dims -dimensions of the sample (meters),
- one, two, or three-element row
- vector
- parameters.npts -number of points in each dimension
- of the sample, one, two, or three-
- element row vector
- parameters.deriv -{'fourier'} requests Fourier diffe-
- rentiation matrices; {'period',n}
