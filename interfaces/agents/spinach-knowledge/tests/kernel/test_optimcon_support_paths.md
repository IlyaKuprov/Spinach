# tests/kernel/test_optimcon_support_paths.m

- Signature: `result=test_optimcon_support_paths()`

## Purpose

Regression-tests small optimal-control support paths: penalties, trapezium-product derivatives, objective collection, and line-search helpers.

## Checks

- A minimal quiet Spinach object supports the low-level calls. `penalty()` is checked for zero value, gradient, and Hessian with no penalty; closed-form norm-square and bounded spillout values, gradients, and Hessians; a derivative norm-square gradient against centred finite differences; and Cartesian amplitude spillout value and gradient against radial references, with its Hessian checked by finite differences.
- `trapdiff()` left- and right-edge derivative matrices are compared with centred finite differences of a small matrix exponential.
- `objeval()` is checked at value, gradient, and Hessian levels against a two-channel objective that subtracts a penalty-like component. The test also checks the corresponding call counters.
- `alpha_conds()` checks monotonic, Armijo, and strong Wolfe curvature acceptance. `cubic_interp()` is checked against a cubic with a halfway maximum. `bracketing()` must accept a short ascent step on a concave quadratic; `sectioning()` must recover its maximum, zero gradient, and successful exit flag.

## Output

- `result` is a regression-test result with explanatory messages.