# kernel/derivatives/fourlap.m

- Signature: `L=fourlap(npoints,extents)`

## Purpose

Returns a Fourier spectral representation of the Laplacian acting on a 3D data array. Syntax: L=fourlap(npoints,extents)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

## Parameters / inputs

- npoints -a three-element vector specifying the number of
- discretization points in each dimension of the
- 3D cube of data that the operator will be acting
- on, ordered as [X Y Z].
- extents -a three-element vector specifying axis extents,
- ordered as [X Y Z].

## Outputs

- L -Fourier spectral Laplacian, a sparse matrix designed to act
- on the vectorization of the 3D data array. The dimensions of
- the data array are assumed to be ordered as [X Y Z].
- Note: periodic boundary conditions.

## Implementation structure

- Returns a Fourier spectral representation of the Laplacian acting
- on a 3D data array. Syntax:
- L=fourlap(npoints,extents)
- npoints - a three-element vector specifying the number of
- discretization points in each dimension of the
- 3D cube of data that the operator will be acting
- on, ordered as [X Y Z].
- extents - a three-element vector specifying axis extents,
- ordered as [X Y Z].
- L -Fourier spectral Laplacian, a sparse matrix designed to act
- on the vectorization of the 3D data array. The dimensions of
- the data array are assumed to be ordered as [X Y Z].
