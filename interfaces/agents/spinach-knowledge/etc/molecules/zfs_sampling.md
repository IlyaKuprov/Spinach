# etc/molecules/zfs_sampling.m

- Signature: `[D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)`

## Purpose

Gadolinium ZFS probability distribution function for DOTA-type ligand complexes in cryogenic water-methanol glasses. The para- meters match those given in Figure 5 of

## Physical / mathematical content

- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Syntax

```matlab
[D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)
```

## Parameters / inputs

- npoints_d -number of Gauss-Legendre quadrature
- points in D
- npoints_e -number of Gauss-Legendre quadrature
- points in E
- tol -tolerance for integration weights
- below which grid points are dropped

## Outputs

- D -a vector of D values at each integration
- grid point
- E -a vector of E values at each integration
- grid point
- W -a vector of weights for each integration
- grid point
- Notes: the function also creates a figure with the distributi-
- ons it has used for D and E parameters.

## Implementation structure

- Gadolinium ZFS probability distribution function for DOTA-type
- ligand complexes in cryogenic water-methanol glasses. The para-
- meters match those given in Figure 5 of
- [D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)
- npoints_d -number of Gauss-Legendre quadrature
- points in D
- npoints_e -number of Gauss-Legendre quadrature
- points in E
- tol -tolerance for integration weights
- below which grid points are dropped
- D -a vector of D values at each integration
- grid point
