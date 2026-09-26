# kernel/grids/gaussleg.m

- Signature: `[x,w]=gaussleg(a,b,n)`

## Purpose

Computes Gauss-Legendre points and weights in [a,b] interval with accuracy order n. Syntax: [x,w]=gaussleg(a,b,n)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

## Parameters / inputs

- a -left edge of the interval
- b -right edge of the interval
- n -accuracy order, the number of points in the
- resulting grid will be n+1.

## Outputs

- x -Gauss-Legendre points
- w -Gauss-Legendre weights

## Implementation structure

- Computes Gauss-Legendre points and weights in [a,b] interval
- with accuracy order n. Syntax:
- [x,w]=gaussleg(a,b,n)
- a -left edge of the interval
- b -right edge of the interval
- n -accuracy order, the number of points in the
- resulting grid will be n+1.
- x -Gauss-Legendre points
- w -Gauss-Legendre weights
- Check consistency
- Initial guess for the nodes in [-1 1]
- Newton-Raphson refinement
