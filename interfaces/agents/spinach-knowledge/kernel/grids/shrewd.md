# kernel/grids/shrewd.m

- Signature: `weights=shrewd(alphas,betas,gammas,max_rank,max_error)`

## Purpose

Computes SHREWD weights for a given two-or three-angle spherical grid. See the paper by Eden and Levitt for details on now the al- gorithm works: http://dx.doi.org/10.1006/jmre.1998.1427 Syntax: weights=shrewd(alphas,betas,gammas,max_rank,max_error)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

## Parameters / inputs

- alphas -alpha Euler angles (ZYZ active) of the
- grid, in radians, set to all-zeros for
- two-angle grids
- betas -beta Euler angles (ZYZ active) of the
- grid, in radians
- gammas -gamma Euler angles (ZYZ active) of the
- grid, in radians
- max_rank -maximum spherical rank to take into consi-
- deration when minimizing residuals
- max_error -maximum residual absolute error per spheri-
- cal function

## Outputs

- weights -a vector of grid weights for each
- [alpha beta gamma] point supplied.

## Implementation structure

- Computes SHREWD weights for a given two-or three-angle spherical
- grid. See the paper by Eden and Levitt for details on now the al-
- gorithm works: http://dx.doi.org/10.1006/jmre.1998.1427 Syntax:
- weights=shrewd(alphas,betas,gammas,max_rank,max_error)
- alphas -alpha Euler angles (ZYZ active) of the
- grid, in radians, set to all-zeros for
- two-angle grids
- betas -beta Euler angles (ZYZ active) of the
- grid, in radians
- gammas -gamma Euler angles (ZYZ active) of the
- max_rank -maximum spherical rank to take into consi-
- deration when minimizing residuals
