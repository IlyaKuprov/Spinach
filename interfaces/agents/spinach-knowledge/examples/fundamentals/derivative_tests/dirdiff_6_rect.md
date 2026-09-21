# examples/fundamentals/derivative_tests/dirdiff_6_rect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/dirdiff_6_rect.m`
- Signature: `dirdiff_6_rect()`
- Total lines: 85

## Purpose

GRAPE Hessian test against finite-differenced gradients.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Implementation structure

- GRAPE Hessian test against finite-differenced gradients.
- Formalisms to test
- Loop over formalisms
- Build the derivative-test system
- Define control parameters
- Spinach housekeeping
- Random guess and finite diff increment
- Call GRAPE and request analytical Hessian
- Leftmost Hessian column
- Rightmost Hessian column
- Middle Hessian column

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dirdiff_test_system()`, `optimcon()`, `eps()`, `grape_xy()`, `squeeze()`, `hess_anl()`, `wave_forw()`, `wave_back()`, `grad_forw()`, `grad_back()`.
