# examples/fundamentals/derivative_tests/dirdiff_5_rect.m

- Signature: `dirdiff_5_rect()`

## Purpose

GRAPE Hessian internal consistency test: Newton against Goodwin algorithm.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

## Implementation structure

- GRAPE Hessian internal consistency test: Newton
- against Goodwin algorithm.
- Formalisms to test
- Loop over formalisms
- Build the derivative-test system
- Define control parameters
- Pick initial guess, phase-modulated GRAPE
- Get Newton Hessian
- Get Goodwin Hessian
- Pick initial guess, XY-modulated GRAPE
- Run the comparisons
