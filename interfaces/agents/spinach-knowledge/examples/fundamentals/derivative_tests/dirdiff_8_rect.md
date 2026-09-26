# examples/fundamentals/derivative_tests/dirdiff_8_rect.m

- Signature: `dirdiff_8_rect()`

## Purpose

GRAPE phase Hessian test against finite-differenced gradients, rectangles integrator.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Implementation structure

- GRAPE phase Hessian test against finite-differenced gradients,
- rectangles integrator.
- Formalisms to test
- Loop over formalisms
- Build the derivative-test system
- Define control parameters
- Set the interval grid
- Spinach housekeeping
- Random phases and finite diff increment
- Call GRAPE and request analytical Hessian
- Leftmost Hessian column
- Rightmost Hessian column
