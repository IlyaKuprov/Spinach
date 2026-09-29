# kernel/derivatives/fdvec.m

Direct source: [kernel/derivatives/fdvec.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fdvec.m)
Spin Dynamics Wiki: [fdvec.m](https://spindynamics.org/wiki/index.php?title=fdvec.m)

## Purpose and interface

dx=fdvec(x,npoints,order) differentiates a numeric row or column vector with a finite-difference stencil and returns derivative estimates in the same shape as x. The function first sizes dx from the input shape, then reshapes a working copy of x to a column for the calculation; this is why row input still produces row output.

- x must be a numeric vector with at least three elements.
- npoints is the positive odd stencil size.
- order is a positive integer smaller than npoints.

## Stencils and coordinate units

For each element near the left edge, fdweights supplies weights evaluated at that element's position among the first npoints samples. The opposite edge uses the reversed weights multiplied by (-1)^order. Between these edges, a centred stencil uses offsets from -(npoints-1)/2 through (npoints-1)/2.

fdvec has no sample-spacing input: the positions passed to fdweights are integer sample indices, so the returned derivative is with respect to that unit-spaced index. It does not impose periodic boundary conditions; the end points are evaluated with sided stencils. For non-unit physical spacing, the caller must account for the corresponding coordinate scaling.

## Input guards and example

The source checks that x is numeric and a vector, that its element count is at least three, that npoints is a positive odd integer, and that order is positive and below npoints. The error text for the length check says “more than three elements,” while the actual predicate rejects only lengths below three; a three-element vector passes that guard. The implementation does not separately check that npoints is no greater than numel(x).

    x=[0 1 4 9 16];
    dx=fdvec(x,3,1);

The call requests a first derivative from five samples using three points per stencil; the output has the same row shape as x.

## Related routine

- [fdweights.m](fdweights.md) computes the coefficients for each stencil.
