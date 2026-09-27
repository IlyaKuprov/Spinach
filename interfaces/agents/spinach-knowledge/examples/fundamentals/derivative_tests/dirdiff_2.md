# examples/fundamentals/derivative_tests/dirdiff_2.m

- Signature: `dirdiff_2()`

## Purpose

Check the analytical left- and right-control derivatives returned by `trapdiff` against central finite differences for the second-order Magnus product quadrature.

## Physical / mathematical content

The test uses a general coherent and non-symmetric dissipative case, represented by separate left and right drift generators and a control operator. The three formalism labels are used to construct Spinach test systems; the derivative check then operates on random matrices.

## Numerical / algorithmic content

The time step is estimated as the mean of the inverse 2-norms of the two drift matrices. The finite-difference increment is `sqrt(eps('double'))`. Analytical directional derivatives are compared with centered differences of matrix exponentials; each difference must be below `10*sqrt(eps('double'))` in 2-norm.

## Implementation structure

- Construct test systems for `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`.
- Generate two random complex `50×50` drift matrices and one random complex control matrix.
- Build the left and right control directions and evaluate both derivatives with `trapdiff`.
- Compare each result with its finite-difference estimate and fail if either check is outside tolerance.
