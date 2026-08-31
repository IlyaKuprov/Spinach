# kernel/operators/lindbladian.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/lindbladian.m`
- Signature: `R=lindbladian(A_left,A_right,rho,rlx_rate)`
- Total lines: 66

## Purpose

Generates a Lindblad superoperator from user-specified left-side and right-side product superoperators and calibrates it using the experi- mental relaxation rate of a user-specified state. Syntax: R=lindbladian(A_left,A_right,rho,rlx_rate)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(A_left,A_right,rho,rlx_rate)`.
- Lines 36-37: Generate a Lindbladian; implemented by `R=A_left*A_right'-(A_left'*A_left+A_right*A_right')/2`.
- Lines 39-40: Check for silly inputs; implemented by `if abs(rho'*R*rho)<1e-10`.
- Lines 44-45: Calibrate the Lindbladian; implemented by `rho=rho/norm(rho,2); R=-rlx_rate*R/(rho'*R*rho)`.

### Control flow inferred from the code

- Line 40: conditional branch on `abs(rho'*R*rho)<1e-10`.

### Key state/data transformations

- Lines 37: computes `R` using `R=A_left*A_right'-(A_left'*A_left+A_right*A_right')/2`.
- Lines 45: computes `rho` using `rho=rho/norm(rho,2); R=-rlx_rate*R/(rho'*R*rho)`.

### Local helper functions

- Line 50: `grumble()` — `function grumble(A_left,A_right,rho,rlx_rate)`. Morality, it could be argued, represents the way that people
  - Representative operation: `if (~isnumeric(A_left))||(~isnumeric(A_right))|| (~isnumeric(rho))||(~isnumeric(rlx_rate))`.
  - Representative operation: `(~isnumeric(rho))||(~isnumeric(rlx_rate))`.

## Parameters / inputs

- A_left -left side product superoperator of the
- interaction that is causing relaxation
- (see operator.m and hamiltonian.m)
- A_right -right side product superoperator of the
- same interaction
- rho -the state vector whose relaxation rate
- is known from the experiment
- rlx_rate -experimental relaxation rate of rho

## Outputs

- R -Lindblad relaxation superoperator indu-
- ced by the interaction A, such that
- <rho|R|rho>/norm(rho,2)^2 = -rlx_rate

## Implementation structure

- Generates a Lindblad superoperator from user-specified left-side and
- right-side product superoperators and calibrates it using the experi-
- mental relaxation rate of a user-specified state. Syntax:
- R=lindbladian(A_left,A_right,rho,rlx_rate)
- A_left -left side product superoperator of the
- interaction that is causing relaxation
- (see operator.m and hamiltonian.m)
- A_right -right side product superoperator of the
- same interaction
- rho -the state vector whose relaxation rate
- is known from the experiment
- rlx_rate -experimental relaxation rate of rho

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
