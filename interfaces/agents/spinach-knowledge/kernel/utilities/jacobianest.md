# kernel/utilities/jacobianest.m

- MATLAB implementation: [kernel/utilities/jacobianest.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/jacobianest.m)

- Signature: `[jac,err] = jacobianest(fun,x0)`
- Source documentation: <https://spindynamics.org/wiki/index.php?title=jacobianest.m>

## Purpose

Estimate the Jacobian of a vector-valued function at a numeric vector or array `x0`, and return an entry-wise estimated error alongside each partial derivative. This is a general numerical-differentiation utility, not a spin-system model.

## Inputs and outputs

- `fun` is a function handle that accepts `x0` and returns a vector-valued result. The implementation checks that it is a function handle.
- `x0` is numeric and may be a vector or array. Its elements are treated as independent scalar coordinates, addressed by linear index; the function receives the original array shape.
- `jac` has one row per element of `fun(x0)` and one column per element of `x0` (the function output is vectorised internally).
- `err` has the same shape as `jac` and contains the estimated error associated with each selected derivative.

## Numerical method

The function evaluates `fun(x0)` once to establish the number of output components. For each input coordinate it forms a geometrically spaced set of 26 perturbations. For a nonzero coordinate, the signed perturbations are `x0(i)*100*(2.0000001).^(0:-1:-25)`; for a zero coordinate, the same relative scale is used without multiplying by zero. Each perturbation gives a centred finite-difference derivative, `(f(x0+h)-f(x0-h))/(2*h)`.

For each output component, `rombextrap` combines successive finite-difference values using the error powers `[2 4]`, cancelling the leading second- and fourth-order error terms, extrapolates towards zero step, and estimates uncertainty from the residual of that extrapolation. The selection step is important: it sorts the *extrapolated derivative values* in ascending numerical order, removes the three smallest and three largest values, and keeps the uncertainty estimates associated with the remaining candidates. It then selects the retained derivative whose estimated error is smallest. It does **not** trim the endpoints of the step-size sequence: the trim acts only after extrapolated derivative candidates and their uncertainties have been formed.

## Shape and edge cases

An input array with `p=numel(x0)` coordinates produces an `n-by-p` Jacobian when `fun(x0)` contains `n` values. If the function result is empty, the routine returns empty `0-by-p` arrays for both outputs. The reported errors are estimates produced by the extrapolation procedure, not an assertion of a rigorous bound.

## Finite differences and error pairing

The code perturbs one input coordinate at a time, evaluates both sides of each centred difference, and then processes each output component independently. `sort` returns indices that are also used to reorder the corresponding error estimates, preserving the derivative/error pairing after the extreme derivative values are discarded.
