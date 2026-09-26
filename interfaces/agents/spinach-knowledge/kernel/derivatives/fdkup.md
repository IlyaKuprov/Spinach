# kernel/derivatives/fdkup.m

- Signature: `K=fdkup(npoints,extents,chi,nstenc)`

## Purpose

Returns a finite difference representation of the Kuprov operator: K[rho]=-(1/3)*Trace(Hessian[rho]*chi) with the number of stencil points in the finite difference approxi- mation specified by user. The resulting operator is a sparse matrix designed to act on the vectorisation of rho. The dimensions of rho are assumed to be ordered as [X Y Z]. For further information, see K=fdkup(npoints,extents,chi,nstenc)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Parameters / inputs

- npoints -a three-element vector specifying the dimensions
- of the 3D cube of data that the operator will be
- acting on, in Angstroms. The dimensions are assu-
- med to be ordered as [X Y Z].
- chi -the electron magnetic susceptibility tensor in
- cubic Angstroms, a symmetric 3x3 matrix.
- extents -a three-element vector specifying axis extents
- in Angstroms. The dimensions are assumed to be
- ordered as [X Y Z].
- nstenc -number of finite-difference stencil points for
- the finite-difference approximation. Periodic
- boundary conditions are used.

## Outputs

- K -a sparse matrix designed to act on the vectori-
- zation of the array. The dimensions are assumed
- to be ordered as [X Y Z].

## Implementation structure

- Returns a finite difference representation of the Kuprov operator:
- K[rho]=-(1/3)*Trace(Hessian[rho]*chi)
- with the number of stencil points in the finite difference approxi-
- mation specified by user. The resulting operator is a sparse matrix
- designed to act on the vectorisation of rho. The dimensions of rho
- are assumed to be ordered as [X Y Z]. For further information, see
- K=fdkup(npoints,extents,chi,nstenc)
- npoints - a three-element vector specifying the dimensions
- of the 3D cube of data that the operator will be
- acting on, in Angstroms. The dimensions are assu-
- med to be ordered as [X Y Z].
- chi - the electron magnetic susceptibility tensor in
